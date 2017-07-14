/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2017 Philip Morris International

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "sectionalFrederix.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    makeFluidThermoTypes(sectionalFrederix, aerosolModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::sectionalFrederix::sectionalFrederix
(
    const fvMesh& mesh
)
:
    aerosolModel(mesh),
    dropletSizeDimension_(dimMass),
    condRateListPtr_(NULL),
    SJDnucPtr_(NULL),
    M0_(0),
    Z0_(0),
    phi_(1.0),
    preparedCoa_(false),
    zStarCoa_(0),
    ijCoa_(0),
    kCoa_(0),
    weightsCoa_(0),
    domainDefect_(0.0),
    zeta_(false),
    mesh_(mesh)
{
    this->read();
    updateSizeDistribution();

    if(sizeDistType_ == NOSIZEDIST)
    {
        FatalErrorIn("Foam::aerosolModels::sectionalFrederix::sectionalFrederix(const fvMesh& mesh)")
            << "This is a sectional method. The size distribution must be other than 'none'." << exit(FatalError);
    }

    M_.setSize(P_);
    V_.setSize(P_);
    J_.setSize(P_);
    S_.setSize(thermo().nSpecies());
    phid_.setSize(P_);

    // Check base field file

    IOobject ioM
    (
        "M",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject ioMP
    (
        word("M."+Foam::name(P_-1)),
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    // Assume that when the N.<P-1> field exists, so do the others

    label l(Foam::log10(Foam::scalar(max(P_-1,1))) + 1);

    if(ioM.headerOk() && !ioMP.headerOk())
    {
        // Create from N and M base fields

        tmp<volScalarField> tM(new volScalarField(ioM, mesh));

        forAll(x_, i)
        {
            Foam::string is(Foam::name(i));
            Foam::string name(std::string(l-is.length(), '0')+is);

            M_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        word("M."+name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tM()
                )
            );
        }
    }
    else if(ioMP.headerOk())
    {
        // Read section N fields directly

        forAll(x_, i)
        {
            Foam::string is(Foam::name(i));
            Foam::string name("M."+std::string(l-is.length(), '0')+is);

            M_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        word(name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
    }
    else
    {
        FatalErrorIn("Foam::aerosolModels::sectionalFrederix::sectionalFrederix(...)")
            << "Could not initialize M fields" << exit(FatalError);
    }

    // Check base field file

    IOobject ioV
    (
        "V",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject ioVP
    (
        word("V."+Foam::name(P_-1)),
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    // Assume that when the V.<P-1> field exists, so do the others

    if(ioV.headerOk() && !ioVP.headerOk())
    {
        // Create from V from V base fields

        tmp<volVectorField> tV(new volVectorField(ioV, mesh));

        forAll(x_, i)
        {
            Foam::string is(Foam::name(i));
            Foam::string name("V."+std::string(l-is.length(), '0')+is);

            V_.set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        word(name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tV()
                )
            );
        }
    }
    else if(ioVP.headerOk())
    {
        // Read section V fields directly

        forAll(x_, i)
        {
            Foam::string is(Foam::name(i));
            Foam::string name("V."+std::string(l-is.length(), '0')+is);

            V_.set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        word(name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
    }
    else
    {
        FatalErrorIn("Foam::aerosolModels::sectionalFrederix::sectionalFrederix(...)")
            << "Could not initialize V fields" << exit(FatalError);
    }

    // Set empty N0 and Z0 fields

    M0_.setSize(P_);

    forAll(x_, i)
    {
        Foam::string is(Foam::name(i));
        Foam::string name("M0."+std::string(l-is.length(), '0')+is);

        M0_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    word(name),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("M0", M_[i].dimensions(), 0.0)
            )
        );
    }

    Z0_.setSize(thermo().nSpecies());

    forAll(thermo().species(), j)
    {
        Z0_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word("Z0." + Foam::name(j)),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("Z0", thermo().Z()[j].dimensions(), 0.0)
            )
        );
    }

    // Set J and S fields to zero fields

    forAll(x_, i)
    {
        Foam::string is(Foam::name(i));
        Foam::string name("J."+std::string(l-is.length(), '0')+is);

        J_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    word(name),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("J", dimDensity*M_[i].dimensions()/dimTime, 0.0)
            )
        );
    }

    forAll(thermo().species(), j)
    {
        S_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word("S." + Foam::name(j)),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("S", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }

    // Set droplet drift flux fields to zero fields

    forAll(x_, i)
    {
        Foam::string is(Foam::name(i));
        Foam::string name("phid."+std::string(l-is.length(), '0')+is);

        phid_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    word(name),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("phid", dimVelocity*dimDensity*dimArea, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::sectionalFrederix::~sectionalFrederix()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::sectionalFrederix::update()
{
    if (doNuc_)
    {
        checkGetNucFieldsOrExit();

        SJDnucPtr_ = NULL;

        SJDnucPtr_ = &getNucFields();
    }

    if (doDrift_)
    {
        checkUpdateDropDriftVelFieldsOrExit();

        updateDropDriftVelFields();

        updateDropletFluxes();
    }

    if (doCond_)
    {
        if (zeta_)
        {
            checkGetCondRateListOrExit();
        }
        else
        {
            checkGetEtaGammaListOrExit();
        }

        condRateListPtr_ = NULL;

        if (zeta_)
        {
            condRateListPtr_ = &getEtaGammaList(x_);
        }
        else
        {
            condRateListPtr_ = &getCondRateList(x_);
        }
    }
}

void Foam::aerosolModels::sectionalFrederix::updateDropletFluxes()
{
    if (doDrift_)
    {
        checkUpdateDropDriftVelFieldsOrExit();

        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        linear<vector> s(mesh_);

        forAll(x_, i)
        {
            phid_[i] = (s.interpolate(thermo().rho() * (V_[i] - U)) & mesh_.Sf());
        }
    }
    else
    {
        forAll(x_, i)
        {
            phid_[i] = dimensionedScalar("phid", dimDensity*dimVelocity*dimArea, 0);
        }
    }
}

void Foam::aerosolModels::sectionalFrederix::fractionalStepInternal()
{
    const dictionary& speciesPhaseChange = thermo().speciesPhaseChange();
    const label& n = thermo().nSpeciesPhaseChange();

    PtrList<volScalarField>& Z = thermo().Z();
    PtrList<volScalarField>& Y = thermo().Y();

    const volScalarField& rho = thermo().rho();

    scalar deltaT = mesh_.time().deltaT().value();

    if (doMonitors_)
    {
        clearScalarMonitors();
    }

    // Condensational growth

    domainDefect_ = 0.0;

    #define DROPCHECK (mag(dz[i]) > 1E-30) && (M0_[i][jCell] > 1.0)
    #define DOMAINCHECK (k[i] >= 0) && (k[i] < P_-1)

    if (doCond_)
    {
        checkGetCondRateListOrExit();

        const PtrList<PtrList<volScalarField> >& Ij = *condRateListPtr_;

        storeM0();
        storeZ0();

        const tmp<volScalarField> tZtot = thermo().Ztot();
        const volScalarField& Ztot = tZtot();

        forAll(mesh_.C(), jCell)
        {
            // Only if we have droplets in this cell

            if (Ztot[jCell] > SMALL)
            {
                List<scalar> zStar(P_, 0.0);
                List<scalar> I(P_, 0.0);

                // Compute total condensation rate and clear current solution

                forAll(x_, i)
                {
                    M_[i][jCell] = 0.0;

                    if (M0_[i][jCell] > 1.0)
                    {
                        forAll(speciesPhaseChange, j)
                        {
                            I[i] += Ij[i][j][jCell];
                        }
                    }
                }

                // Compute new size

                if (zeta_)
                {
                    List<scalar> zetaStar(P_, 0.0);

                    forAll(x_, i)
                    {
                        zetaStar[i] = Ij[i][n][jCell] + I[i]*deltaT;
                    }

                    // Compute real z. PsiInv() gives -1 if zetaStar < 0.

                    zStar = psiInv(zetaStar, jCell);
                }
                else
                {
                    forAll(x_, i)
                    {
                        zStar[i] = x_[i] + I[i]*deltaT;
                    }
                }

                // Let droplets not grow beyond the sectional domain

                forAll(x_, i)
                {
                    if (zStar[i] > x_[P_-1])
                    {
                        zStar[i] = x_[i];

                        domainDefect_ += M0_[i][jCell] * x_[i] * rho[jCell];
                    }
                }

                // Update species mass fractions

                List<scalar> dz = zStar - x_;

                // Find the index of the section which is directly left to zStar

                List<label> k = xLowerIndex(zStar);

                forAll(x_, i)
                {
                    if (DROPCHECK)
                    {
                        if (DOMAINCHECK)
                        {
                            // Evaporation/condensation within sections

                            forAll(speciesPhaseChange, j)
                            {
                                scalar W = Ij[i][j][jCell] / stabilise(I[i], 1E-99);
                                scalar S = dz[i] * W * M0_[i][jCell];

                                Z[j][jCell] += S;
                                Y[j][jCell] -= S;

                                if (doMonitors_)
                                {
                                    (*scalarMonitorPtrs_[word("S." + Foam::name(j))])[jCell] += S/deltaT;
                                }
                            }
                        }
                        else if (zStar[i] < x_[0])
                        {
                            // Complete evaporation

                            scalar sumZ = 0.0;

                            forAll(speciesPhaseChange, j)
                            {
                                sumZ += Z0_[j][jCell];
                            }

                            forAll(speciesPhaseChange, j)
                            {
                                scalar S = x_[i] * M0_[i][jCell] * Z0_[j][jCell] / stabilise(sumZ, 1E-99);

                                Z[j][jCell] -= S;
                                Y[j][jCell] += S;

                                if (doMonitors_)
                                {
                                    (*scalarMonitorPtrs_[word("S." + Foam::name(j))])[jCell] += S/deltaT;
                                }
                            }
                        }
                    }
                }

                // Compute new solution

                switch(distMethod_)
                {
                    case TWOMOMENT:
                    {
                        forAll(x_, i)
                        {
                            if (DOMAINCHECK)
                            {
                                if (DROPCHECK && DOMAINCHECK)
                                {
                                    twoMoment(k[i], jCell, zStar[i], M0_[i][jCell]);
                                }
                                else
                                {
                                    M_[i][jCell] += M0_[i][jCell];
                                }
                            }
                        }
                    }
                    break;

                    case FOURMOMENT:
                    {
                        forAll(x_, i)
                        {
                            if (DOMAINCHECK)
                            {
                                if (DROPCHECK)
                                {
                                    fourMoment(k[i], jCell, zStar[i], M0_[i][jCell]);
                                }
                                else
                                {
                                    M_[i][jCell] += M0_[i][jCell];
                                }
                            }
                        }
                    }
                    break;

                    case HYBRID:
                    {
                        // Compute four-moment hybrid solution

                        forAll(x_, i)
                        {
                            if (DOMAINCHECK)
                            {
                                if (DROPCHECK)
                                {
                                    fourMoment(k[i], jCell, zStar[i], M0_[i][jCell], phi_);
                                    twoMoment(k[i], jCell, zStar[i], M0_[i][jCell], 1.0-phi_);
                                }
                                else
                                {
                                    M_[i][jCell] += M0_[i][jCell];
                                }
                            }
                        }

                        // Revert negative values to two-moment

                        List<label> revert(P_, false);

                        forAll(x_, i)
                        {
                            if (M_[i][jCell] < 0.0)
                            {
                                forAll(x_, j)
                                {
                                    if (k[j] >= i-2 && k[j] <= i+1)
                                    {
                                        revert[j] = true;
                                    }
                                }
                            }
                        }

                        forAll(x_, i)
                        {
                            if (revert[i] && DROPCHECK && DOMAINCHECK)
                            {
                                fourMoment(k[i], jCell, zStar[i], M0_[i][jCell], -phi_);
                                twoMoment(k[i], jCell, zStar[i], M0_[i][jCell], phi_);
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    // Nucleation

    if (doNuc_)
    {
        checkGetNucFieldsOrExit();

        const PtrList<volScalarField>& SJDnuc = *SJDnucPtr_;

        forAll(mesh_.C(), jCell)
        {
            if (SJDnuc[n][jCell] > 0)
            {
                scalar z = SJDnuc[n+1][jCell];

                if (doMonitors_)
                {
                    (*scalarMonitorPtrs_["znuc"])[jCell] = SJDnuc[n+1][jCell];
                }

                if (z < x_[0])
                {
                    M_[0][jCell] += deltaT*SJDnuc[n][jCell]/rho[jCell] * z/x_[0];

                    forAll(speciesPhaseChange, j)
                    {
                        Z[j][jCell] += deltaT*SJDnuc[j][jCell]/rho[jCell];
                        Y[j][jCell] -= deltaT*SJDnuc[j][jCell]/rho[jCell];
                    }

                    if (doMonitors_)
                    {
                        (*scalarMonitorPtrs_["Jnuc"])[jCell] = SJDnuc[n][jCell];
                    }
                }
                else if (z >= x_[P_-1])
                {
                    FatalErrorIn("Foam::aerosolModels::sectionalFrederix::fractionalStep()")
                        << "Nucleation is occuring outside the size domain" << exit(FatalError);
                }
                else
                {
                    label i = xLowerIndex(z);
                    twoMoment(i, jCell, z, deltaT*SJDnuc[n][jCell]/rho[jCell]);

                    forAll(speciesPhaseChange, j)
                    {
                        Z[j][jCell] += deltaT*SJDnuc[j][jCell]/rho[jCell];
                        Y[j][jCell] -= deltaT*SJDnuc[j][jCell]/rho[jCell];
                    }

                    if (doMonitors_)
                    {
                        (*scalarMonitorPtrs_["Jnuc"])[jCell] = SJDnuc[n][jCell];
                    }
                }
            }
        }
    }

    // Update boundaries

    if (doNuc_ || doCond_)
    {
        forAll(x_, i)
        {
            M_[i].correctBoundaryConditions();
        }

        forAll(thermo().speciesPhaseChange(), j)
        {
            Y[j].correctBoundaryConditions();
            Z[j].correctBoundaryConditions();
        }
    }
}

void Foam::aerosolModels::sectionalFrederix::fractionalStepExternal()
{
    if (doMonitors_)
    {
        *scalarMonitorPtrs_["dcm"] = dcm();
        *scalarMonitorPtrs_["dmm"] = dmm();
    }

    if (doCoa_)
    {
        checkGetCoaRateFieldOrExit();

        if(!preparedCoa_)
        {
            prepareCoa();
        }

        scalar deltaT = mesh_.time().deltaT().value();

        const volScalarField& p1 = thermo().p1();
        const dimensionedScalar p0 = thermo().p0();
        const volScalarField& T = thermo().T();
        const volScalarField& rho = thermo().rho();

        const dictionary speciesPhaseChange = thermo().speciesPhaseChange();

        const dictionary& species = thermo().species();
        const label& nSpecies = thermo().nSpecies();

        const List<Foam::scalar>& M = thermo().M();

        const PtrList<Foam::volScalarField>& Y = thermo().Y();

        const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

        const scalar pi = constant::mathematical::pi;
        const scalar kB = 1.3806488E-23;
        const scalar N_A = 6.0221413E+23;

        // Molecular mass

        List<Foam::scalar> m(nSpecies);

        forAll(species, i)
        {
            m[i] = 0.001*M[i]/N_A;
        }

        // Coalescence

        storeM0();

        List<scalar> v(P_, 0.0);

        // Get liquid density

        tmp<volScalarField> tRhoLiquid = thermo().rhoLiquid();
        volScalarField& rhoLiquid = tRhoLiquid();

        forAll(mesh_.C(), jCell)
        {
            if (rhoLiquid[jCell] > 0.0)
            {
                // We have droplets. Compute droplet volumes

                forAll(x_, i)
                {
                    v[i] = x_[i]/rhoLiquid[jCell];
                }

                // Compute mean free path

                scalar sumY = 0.0;
                scalar sumYm = 0.0;

                forAll(species, j)
                {
                    sumY += Y[j][jCell];
                    sumYm += Y[j][jCell]/m[j];
                }

                scalar mg = sumY/sumYm;

                scalar lambda = sqrt(8.0*kB*T[jCell]/pi/mg)*(4.0*muEff[jCell]/5.0/(p1[jCell]+p0.value()));

                // Loop over all possible combinations

                forAll(kCoa_, l)
                {
                    scalar i = ijCoa_[l][0];
                    scalar j = ijCoa_[l][1];

                    if (M0_[i][jCell] > SMALL && M0_[j][jCell] > SMALL)
                    {
                        scalar a = pow(v[i], 1.0/3.0) + pow(v[j], 1.0/3.0);

                        scalar d = pow(6.0/(8.0*pi), 1.0/3.0) * a;

                        scalar Kn = lambda/d;

                        scalar fij = M0_[i][jCell] * M0_[j][jCell] * rho[jCell] *
                            getCoaRateCell(v[i], v[j], jCell, Kn) * deltaT;

                        scalar k = kCoa_[l];

                        M_[i][jCell] -= fij;
                        M_[j][jCell] -= fij;

                        M_[k][jCell] += weightsCoa_[l][0] * fij;
                        M_[k+1][jCell] += weightsCoa_[l][1] * fij;
                    }
                }
            }
        }

        // Update boundaries

        forAll(x_, i)
        {
            M_[i].correctBoundaryConditions();
        }

        if (doMonitors_)
        {
            // Approximate coa rate

            forAll(x_, i)
            {
                *scalarMonitorPtrs_["Jcoa"] += rho * (M_[i]-M0_[i])/mesh_.time().deltaT();
            }
        }
    }
}

void Foam::aerosolModels::sectionalFrederix::checkConsistency()
{
    const PtrList<volScalarField>& Z = thermo().Z();
    const volScalarField& rho = thermo().rho();

    scalar MNsum = 0.0;
    scalar MZsum = 0.0;

    forAll(thermo().species(), j)
    {
        MZsum +=
            sum
            (
                rho.internalField() * Z[j].internalField() * mesh_.V().field()
            );
    }

    forAll(x_, i)
    {
        MNsum +=
            sum
            (
                x_[i] * M_[i].internalField() * rho.internalField() * mesh_.V().field()
            );
    }

    reduce(MNsum, sumOp<scalar>());
    reduce(MZsum, sumOp<scalar>());

    reduce(domainDefect_, sumOp<scalar>());

    Info << "Total liquid mass: " << MZsum << " kg (Zj), " << MNsum << " kg (Ni), difference = "
         << mag(MNsum-MZsum) << ", domain defect = "
         << (domainDefect_/stabilise(MZsum, 1E-99)*100) << "%" << endl;
}

void Foam::aerosolModels::sectionalFrederix::correctSizeDistribution()
{
    if (doCorrSizeDist_)
    {
        const PtrList<volScalarField>& Z = thermo().Z();

        forAll(mesh_.C(), iCell)
        {
            scalar MN = 0.0;

            forAll(x_, i)
            {
                MN += M_[i][iCell] * x_[i];
            }

            scalar MZ = 0.0;

            forAll(thermo().species(), j)
            {
                MZ += Z[j][iCell];
            }

            // Adjust M to match mass

            forAll(x_, i)
            {
                M_[i][iCell] *= MZ/stabilise(MN, 1E-99);
            }
        }

        forAll(x_, i)
        {
            M_[i].correctBoundaryConditions();
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::sectionalFrederix::dcm()
{
    const scalar pi = constant::mathematical::pi;

    tmp<volScalarField> tdcm
    (
        new volScalarField
        (
            IOobject
            (
                "dcm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dcm", dimLength, 0.0)
        )
    );

    volScalarField& dcm = tdcm();

    volScalarField temp1
    (
        IOobject
        (
            "temp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("temp", dimLength/dimMass, 0.0)
    );

    volScalarField temp2
    (
        IOobject
        (
            "temp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("temp", dimless/dimMass, 0.0)
    );

    tmp<volScalarField> tRhol = thermo().rhoLiquid();
    volScalarField& rhol = tRhol();

    thermo().limitLiquidDensity(rhol);

    dimensionedScalar unity("small", dimless/dimMass, 1.0);

    forAll(x_, i)
    {
        temp1 += pow(dimensionedScalar("z", dimMass, x_[i]), 1.0/3.0) * M_[i]
               * pow(6.0/pi, 1.0/3.0)
               / pow(rhol, 1.0/3.0);

        temp2 += M_[i];
    }

    dcm == temp1/stabilise(temp2, unity);

    dcm.max(dimensionedScalar("zero", dimLength, 0.0));

    return tdcm;
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::sectionalFrederix::dmm()
{
    const scalar pi = constant::mathematical::pi;

    tmp<volScalarField> tdmm
    (
        new volScalarField
        (
            IOobject
            (
                "dmm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dmm", dimLength, 0.0)
        )
    );

    volScalarField& dmm = tdmm();

    volScalarField temp1
    (
        IOobject
        (
            "temp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("temp", dimLength, 0.0)
    );

    volScalarField temp2
    (
        IOobject
        (
            "temp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("temp", dimless, 0.0)
    );

    tmp<volScalarField> tRhol = thermo().rhoLiquid();
    volScalarField& rhol = tRhol();

    thermo().limitLiquidDensity(rhol);

    dimensionedScalar small("small", dimless, SMALL);

    forAll(x_, i)
    {
        temp1 += pow(dimensionedScalar("z", dimMass, x_[i]), 4.0/3.0) * M_[i]
               * pow(6.0/pi, 1.0/3.0)
               / pow(rhol, 1.0/3.0);

        temp2 += dimensionedScalar("z", dimMass, x_[i]) * M_[i];
    }

    dmm == temp1/stabilise(temp2, small);

    dmm.max(dimensionedScalar("zero", dimLength, 0.0));

    return tdmm;
}

void Foam::aerosolModels::sectionalFrederix::storeM0()
{
    forAll(x_, i)
    {
        M0_[i] == M_[i];
    }
}

void Foam::aerosolModels::sectionalFrederix::storeZ0()
{
    forAll(thermo().species(), j)
    {
        Z0_[j] == thermo().Z()[j];
    }
}

void Foam::aerosolModels::sectionalFrederix::prepareCoa()
{
    // Pre-compute coalescence lists

    forAll(x_, i)
    {
        for(label j = i; j < P_; j++)
        {
            scalar zStar = x_[i] + x_[j];

            if (zStar < x_[P_-1])
            {
                zStarCoa_.append(zStar);

                ijCoa_.append(List<label>(2, 0));
                ijCoa_[ijCoa_.size()-1][0] = i;
                ijCoa_[ijCoa_.size()-1][1] = j;

                scalar k = xLowerIndex(zStar);

                kCoa_.append(k);

                weightsCoa_.append(weights(x_[k], x_[k+1], zStar));
            }
        }
    }

    preparedCoa_ = true;
}

bool Foam::aerosolModels::sectionalFrederix::read()
{
    if (aerosolModel::read())
    {
        distMethod_ = distMethodNames.read(params_.lookup("distMethod"));

        if (distMethod_ == HYBRID)
        {
            params_.lookup("phi") >> phi_;

            if ((phi_ < 0.0) | (phi_ > 1.0))
            {
                FatalErrorIn("Foam::aerosolModels::sectionalFrederix::read()")
                    << "The value for phi must lie in between or be equal to 0 and 1." << exit(FatalError);
            }
        }

        params_.lookup("solveInZeta") >> zeta_;

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
