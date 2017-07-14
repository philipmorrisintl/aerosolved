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

#include "twoMomentLogNormalFrederix.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    makeFluidThermoTypes(twoMomentLogNormalFrederix, aerosolModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormalFrederix::twoMomentLogNormalFrederix
(
    const fvMesh& mesh
)
:
    aerosolModel(mesh),
    dropletSizeDimension_(dimless),
    tDistData_(2),
    mesh_(mesh)
{
    this->read();
    updateSizeDistribution();

    if(sizeDistType_ != NOSIZEDIST)
    {
        FatalErrorIn("Foam::aerosolModels::twoMomentLogNormalFrederix::twoMomentLogNormalFrederix(const fvMesh& mesh)")
            << "This is a moment method. The size distribution type must be 'none'." << exit(FatalError);
    }

    M_.setSize(1);
    V_.setSize(1);
    J_.setSize(1);
    S_.setSize(thermo().nSpecies());
    phid_.setSize(1);

    M_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "M",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    V_.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "V",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    J_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "J",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("J", dimless/dimVolume/dimTime, 0.0)
        )
    );

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

    tDistData_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "CMD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("CMD", dimLength, 0.0),
            "zeroGradient"
        )
    );

    tDistData_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "Q",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Q", dimless/dimVolume, 0.0)
        )
    );

    phid_.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                "phid",
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormalFrederix::~twoMomentLogNormalFrederix()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::twoMomentLogNormalFrederix::updateDistData()
{
    const volScalarField& rho = thermo().rho();

    PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo().getProperty("rho", fluidThermo::LIQUID);

    const scalar pi = constant::mathematical::pi;

    volScalarField& CMD = tDistData_[0];
    volScalarField& Q = tDistData_[1];

    CMD *= 0.0;
    Q *= 0.0;

    forAll(mesh_.C(), iCell)
    {
        // Compute liquid volume fraction W

        scalar W = 0.0;

        forAll(thermo().species(), j)
        {
            scalar rho_l = dataEntriesListRho_l[j].value(min(thermo().T()[iCell],thermo().Tc()[j]));
            W += thermo().Z()[j][iCell]/stabilise(rho_l, SMALL);
        }

        W *= thermo().rho()[iCell];

        if (W > 0.0 && M_[0][iCell] > 1.0)
        {
            // Initial guess for zero Dcrit

            CMD[iCell] = pow(6.0*W/(pi*M_[0][iCell]*rho[iCell]), 1.0/3.0) * exp(-3.0/2.0*sqr(log(W_)));

            // If the initial guess CMD is smaller than Dcrit_, use zero Dcrit
            // and the initial guess as solution

            if (CMD[iCell] <= Dcrit_)
            {
                Q[iCell] = M_[0][iCell] * rho[iCell];
            }
            else
            {
                scalar alpha = pi*M_[0][iCell] * rho[iCell] / (6.0 * W) * exp(9.0/2.0 * sqr(log(W_)));
                scalar beta = 1/(sqrt(2.0)*log(W_));
                scalar gamma = log(Dcrit_);
                scalar delta = 3.0*sqr(log(W_));

                scalar f = 0.0;
                scalar dfdx = 0.0;
                scalar g = 0.0;
                scalar h = 0.0;

                label Niter;

                for (Niter=1; Niter<=maxIter_; Niter++)
                {
                    g = beta*(gamma - log(CMD[iCell]) - delta);
                    h = beta*(gamma - log(CMD[iCell]));
                    f = pow(CMD[iCell], 3.0) * alpha * erfc(g)/erfc(h) - 1.0;
                    dfdx = alpha * 3 * sqr(CMD[iCell]) * erfc(g) / erfc(h)
                         + alpha * pow(CMD[iCell], 2.0) / erfc(h) * 2.0/sqrt(pi) * exp(-sqr(g))
                         - alpha * pow(CMD[iCell], 2.0) * erfc(g) / sqr(erfc(h)) * 2.0/sqrt(pi) * exp(-sqr(h));

                    CMD[iCell] -= f/dfdx;

                    if (mag(f/dfdx)/mag(CMD[iCell]) < TOL_) break;
                }

                Q[iCell] = M_[0][iCell] * rho[iCell] * 2.0 / (1.0 - erf(log(Dcrit_/CMD[iCell])/(sqrt(2.0) * log(W_))));
            }
        }
    }
    CMD.correctBoundaryConditions();
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::update()
{
    const volScalarField& rho = thermo().rho();

    PtrList<volScalarField>& SJDnuc = getNucFields();

    const scalar pi = constant::mathematical::pi;
    const scalar kB = 1.3806488E-23;
    const scalar N_A = 6.0221413E+23;

    const List<scalar> MM = thermo().M();

    PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo().getProperty("rho", fluidThermo::LIQUID);

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

    const label n = thermo().nSpeciesPhaseChange();

    const volScalarField& T = thermo().T();
    const volScalarField& p1 = thermo().p1();
    const dimensionedScalar p0 = thermo().p0();

    if (doMonitors_)
    {
        clearScalarMonitors();
    }

    // Molecular mass

    List<scalar> m(thermo().nSpecies());

    forAll(thermo().species(), j)
    {
        m[j] = 0.001*MM[j]/N_A;
    }

    // Set to zero

    J_[0] *= 0.0;

    forAll(thermo().speciesPhaseChange(), j)
    {
        S_[j] *= 0.0;
    }

    updateDistData();

    // Update drift

    updateDropletFluxes();

    // Condensation and evaporation mass transport

    volScalarField& CMD = tDistData_[0];
    volScalarField& Q = tDistData_[1];

    List<scalar> x(2*Mint_+1, 0.0);
    List<scalar> d(2*Mint_+1, 0.0);
    List<scalar> z(2*Mint_+1, 0.0);

    forAll(mesh_.C(), iCell)
    {
        // Only do something if the upper boundary of the integration domain
        // (four times the geometric standard deviation) is larger than the
        // critical diameter

        if (doCond_ && CMD[iCell]*pow(W_, 4.0) > Dcrit_)
        {
            scalar dx = (log(CMD[iCell]) + 4.0*log(W_) - log(Dcrit_))/Mint_;

            x[0] = log(Dcrit_);

            for(label i = 1; i <= Mint_; i++)
            {
                x[i*2] = log(Dcrit_) + dx * i;
                x[i*2-1] = (x[i*2]+x[i*2-2])/2.0;
            }

            forAll(x, i)
            {
                d[i] = exp(x[i]);
            }

            scalar sumZ = 0.0;
            scalar sumZrho = 0.0;

            forAll(thermo().speciesPhaseChange(), j)
            {
                sumZ += max(thermo().Z()[j][iCell],0.0);
            }

            forAll(thermo().speciesPhaseChange(), j)
            {
                sumZrho += max(thermo().Z()[j][iCell],0.0)/dataEntriesListRho_l[j].value(min(thermo().T()[iCell],thermo().Tc()[j]));
            }

            scalar rhoMean = sumZ/stabilise(sumZrho, 1E-99);

            // Compute droplet mass

            forAll(z, i)
            {
                z[i] = pow(d[i],3.0) * rhoMean * pi / 6.0;
            }

            // Condensation/evaporation rate in kg/s

            List<List<scalar> > I = getCondRateListCell(z, iCell);

            scalar sumIcrit = 0.0;

            forAll(thermo().speciesPhaseChange(), j)
            {
                S_[j][iCell] = 0.0;

                List<scalar> f(d.size());

                forAll(x, i)
                {
                    f[i] = I[i][j] * Q[iCell] / (sqrt(2.0*pi) * log(W_))
                         * exp(-sqr(x[i] - log(CMD[iCell]))/(2.0*sqr(log(W_))));
                }

                // Simpson's rule

                for(label i = 1; i <= Mint_; i++)
                {
                    S_[j][iCell] += dx/6.0 * (f[i*2-2] + 4.0*f[i*2-1] + f[i*2]);
                }

                sumIcrit += I[0][j];
            }

            if (sumIcrit < 0.0)
            {
                // Complete evaporation

                J_[0][iCell] = sumIcrit*2.0 / (rhoMean * sqr(Dcrit_) * pi)
                             * Q[iCell] / (sqrt(2.0*pi) * Dcrit_ * log(W_))
                             * exp(-sqr(log(Dcrit_) - log(CMD[iCell]))/(2.0*sqr(log(W_))));
            }
        }

        // Nucleation mass transport and rate (only if the nuc diameter > Dcrit)

        if (doNuc_ && SJDnuc[n+1][iCell] > 0.0)
        {
            scalar sumZ = 0.0;
            scalar sumZrho = 0.0;

            forAll(thermo().speciesPhaseChange(), j)
            {
                sumZ += SJDnuc[j][iCell];
            }

            forAll(thermo().speciesPhaseChange(), j)
            {
                sumZrho += SJDnuc[j][iCell]/dataEntriesListRho_l[j].value(min(thermo().T()[iCell],thermo().Tc()[j]));
            }

            scalar rhoMean = sumZ/stabilise(sumZrho, 1E-99);

            scalar dnuc = pow(SJDnuc[n+1][iCell] / stabilise(rhoMean, SMALL) * 6.0/pi, 1.0/3.0);

            if (doMonitors_)
            {
                (*scalarMonitorPtrs_["znuc"])[iCell] = SJDnuc[n+1][iCell];
            }

            if (dnuc > Dcrit_)
            {
                forAll(thermo().speciesPhaseChange(), j)
                {
                    S_[j][iCell] += SJDnuc[j][iCell];
                }

                J_[0][iCell] += SJDnuc[n][iCell];

                if (doMonitors_)
                {
                    (*scalarMonitorPtrs_["Jnuc"])[iCell] += J_[0][iCell];
                }
            }
        }

        // Coalescence rate (only if CMD is larger than Dcrit, and if we have
        // droplets)

        if (doCoa_ && CMD[iCell] > Dcrit_ && M_[0][iCell] > SMALL)
        {
            // Compute volume of droplet

            scalar v = pi/6.0 * pow(CMD[iCell], 3.0);

            // Compute mean free path and Knudsen number

            scalar sumY = 0.0;
            scalar sumYm = 0.0;

            forAll(thermo().species(), j)
            {
                sumY += max(thermo().Y()[j][iCell],0.0);
                sumYm += max(thermo().Y()[j][iCell],0.0)/m[j];
            }

            scalar mg = sumY/sumYm;

            scalar lambda = sqrt(8.0*kB*T[iCell]/pi/mg)*(4.0*muEff[iCell]/5.0/(p1[iCell]+p0.value()));

            scalar Kn = lambda*2.0/CMD[iCell];

            // Get coalesence rate

            scalar beta = getCoaRateCell(v, v, iCell, Kn);

            J_[0][iCell] -= beta * sqr(max(M_[0][iCell],0.0) * rho[iCell]);

            if (doMonitors_)
            {
                (*scalarMonitorPtrs_["Jcoa"])[iCell] -= beta * sqr(max(M_[0][iCell],0.0) * rho[iCell]);
            }
        }
    }

    if (doMonitors_)
    {
        forAll(thermo().speciesPhaseChange(), j)
        {
            (*scalarMonitorPtrs_[word("S." + Foam::name(j))]) == S_[j];
        }
    }

    // Update enthalpy of evaporation

    updateHvapS();
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::updateDropletFluxes()
{
    if (doDrift_)
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        const volScalarField dMeanA(tDistData_[0] * exp(sqr(log(W_))));

        linear<vector> s(mesh_);

        checkUpdateDropDriftVelFieldsOrExit();

        updateDropDriftVelField(dMeanA);

        phid_[0] = (s.interpolate(thermo().rho() * (V_[0] - U)) & mesh_.Sf());
    }
    else
    {
        phid_[0] = dimensionedScalar("phid", dimDensity*dimVelocity*dimArea, 0);
    }
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::fractionalStepInternal()
{
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::fractionalStepExternal()
{
    if (doMonitors_)
    {
        *scalarMonitorPtrs_["dcm"] = dcm();
        *scalarMonitorPtrs_["dmm"] = dmm();
    }
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::checkConsistency()
{
}

void Foam::aerosolModels::twoMomentLogNormalFrederix::correctSizeDistribution()
{
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::twoMomentLogNormalFrederix::dcm()
{
    updateDistData();

    tmp<volScalarField> tdcm
    (
        new volScalarField
        (
            "dcm",
            tDistData_[0] * exp(0.5*sqr(log(W_)))
        )
    );

    tdcm().max(dimensionedScalar("zero", tdcm().dimensions(), 0.0));

    return tdcm;
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::twoMomentLogNormalFrederix::dmm()
{
    updateDistData();

    tmp<volScalarField> tdmm
    (
        new volScalarField
        (
            "dmm",
            tDistData_[0] * exp(3.5*sqr(log(W_)))
        )
    );

    tdmm().max(dimensionedScalar("zero", tdmm().dimensions(), 0.0));

    return tdmm;
}

bool Foam::aerosolModels::twoMomentLogNormalFrederix::read()
{
    if (aerosolModel::read())
    {
        params_.lookup("W") >> W_;
        params_.lookup("M") >> Mint_;
        params_.lookup("Dcrit") >> Dcrit_;
        params_.lookup("maxIter") >> maxIter_;
        params_.lookup("TOL") >> TOL_;

        thermo().readProperty("rho", fluidThermo::LIQUID, thermo().species());

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
