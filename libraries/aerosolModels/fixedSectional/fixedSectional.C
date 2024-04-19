/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2019 Philip Morris International

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

#include "addToRunTimeSelectionTable.H"
#include "fixedSectional.H"
#include "fv.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    defineTypeNameAndDebug(fixedSectional, 0);
    addToRunTimeSelectionTable(aerosolModel, fixedSectional, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::aerosolModels::fixedSectional::updateDrift()
{
    phiInertial_ *= 0.0;
    phiBrownian_ *= 0.0;
    DDisp_ *= 0.0;

    const volScalarField alpha(max(system_->alpha(), residualAlpha_));

    if (this->drift().inertial().type() != "none")
    {
        const volScalarField& rho = this->rho();

        forAll(system_->distribution(), i)
        {
            const word sectionName(system_->distribution()[i].sectionName());

            surfaceScalarField& phiInertiali =
                system_->distribution()[i].phiInertial();

            const volScalarField& Mi = system_->distribution()[i].M();

            const volScalarField d(system_->d(i));

            phiInertiali = fvc::flux
            (
                this->drift().inertial().V(d, sectionName)*rho
            );

            phiInertial_ +=
                system_->distribution()[i].xd()
              * phiInertiali
              * linearInterpolate(Mi);
        }

        phiInertial_ /= linearInterpolate(alpha);
    }
    else
    {
        forAll(system_->distribution(), i)
        {
            system_->distribution()[i].phiInertial() *= 0.0;
        }
    }

    if (this->drift().Brownian().type() != "none")
    {
        const volScalarField& rho = this->rho();

        tmp<volVectorField> tUDiff
        (
            new volVectorField
            (
                IOobject
                (
                    "UDiff",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("UDiff", dimVelocity, vector::zero)
            )
        );

        volVectorField& UDiff = tUDiff.ref();

        forAll(system_->distribution(), i)
        {
            surfaceScalarField& phiBrowniani =
                system_->distribution()[i].phiBrownian();

            const volScalarField& Mi = system_->distribution()[i].M();

            const volScalarField d(system_->d(i));

            volScalarField& DDispi = system_->distribution()[i].D();

            DDispi = this->drift().Brownian().D(d);

            phiBrowniani = -fvc::snGrad(rho*DDispi)*mesh_.magSf();

            UDiff +=
                system_->distribution()[i].xd()
              * DDispi
              * fvc::grad(Mi);

            phiBrownian_ +=
                system_->distribution()[i].xd()
              * phiBrowniani
              * linearInterpolate(Mi);
        }

        const scalar delta(gAverage(Foam::pow(mesh_.V().field(), 1.0/3.0)));

        const dimensionedScalar smallGradAlpha
        (
            "grad(alpha)",
            alpha.dimensions()/dimLength,
            residualAlpha_.value()/delta
        );

        const volScalarField magGradAlpha
        (
            max(mag(fvc::grad(alpha)), smallGradAlpha)
        );

        DDisp_ = mag(UDiff)/magGradAlpha;

        phiBrownian_ /= linearInterpolate(alpha);
    }
    else
    {
        forAll(system_->distribution(), i)
        {
            system_->distribution()[i].D() *= 0.0;
        }
    }
}

void Foam::aerosolModels::fixedSectional::solveSpatial()
{
    Info<<"fixedSectional: solving spatial step" << endl;

    const volScalarField& rho = this->rho();

    const surfaceScalarField& phi = this->phi();

    // Compute the relative and corrective sectional fluxes

    PtrList<surfaceScalarField> phiRelM(system_->distribution().size());

    surfaceScalarField phiCorr(phi*0.0);

    forAll(system_->distribution(), i)
    {
        const volScalarField& Mi = system_->distribution()[i].M();

        const surfaceScalarField& phiInertiali =
            system_->distribution()[i].phiInertial();

        const surfaceScalarField& phiBrowniani =
            system_->distribution()[i].phiBrownian();

        const surfaceScalarField phiReli
        (
            phiInertiali
          + phiBrowniani
          - phiInertial_
          - phiBrownian_
        );

        tmp<fv::convectionScheme<scalar>> convPhiReli
        (
            fv::convectionScheme<Foam::scalar>::New
            (
                mesh_,
                phiReli,
                mesh_.divScheme("div(mvConv)")
            )
        );

        phiRelM.set
        (
            i,
            new surfaceScalarField(convPhiReli->flux(phiReli, Mi))
        );

        phiCorr -= system_->distribution()[i].xd()*phiRelM[i];
    }

    phiCorr /= max(linearInterpolate(system_->alpha()), SMALL);

    // Solve the system of equations

    forAll(system_->distribution(), i)
    {
        volScalarField& Mi = system_->distribution()[i].M();

        const surfaceScalarField phiCorrM(phiCorr*linearInterpolate(Mi));

        const volScalarField D
        (
            rho*system_->distribution()[i].D()
          + turbulence().mut()
        );

        // TODO: is creation here really a good idea?
        fv::options& fvOptions(fv::options::New(mesh_));

	// Number concentration Equation
        fvScalarMatrix MEqn
        (
            fvm::ddt(rho, Mi)
          + mvPhi_->fvmDiv(phi, Mi)
          + mvPhiInertial_->fvmDiv(phiInertial_, Mi)
          + mvPhiBrownian_->fvmDiv(phiBrownian_, Mi)
          - mvPhiDrift_->fvmDiv(phiDrift_, Mi)
          + fvc::div(phiRelM[i]+phiCorrM)
          ==
            fvm::laplacian(D, Mi,"laplacian(D,Mi)")
          + fvOptions(rho, Mi)
        );

        MEqn.relax();

        fvOptions.constrain(MEqn);

        MEqn.solve(mesh_.solver("M"));

        fvOptions.correct(Mi);

        Mi.max(0.0);

        phiEff_[i] = MEqn.flux() + phiRelM[i] + phiCorrM;
    }

    system_->rescale();
}

void Foam::aerosolModels::fixedSectional::solveInternal()
{
    clearRates();

    if
    (
        nucleation_->modelType() == "none"
     && condensation_->modelType() == "none"
     && coalescence_->modelType() == "none"
    )
    {
        return;
    }

    Info<<"fixedSectional: solving internal step" << endl;

    //const scalar pi = constant::mathematical::pi;

    const speciesTable& activeSpecies = thermo_.activeSpecies();
    const speciesTable& contSpecies = thermo_.contSpecies();

    const scalarField& p = thermo_.p().field();
    const scalarField& T = thermo_.T().field();
    const scalarField& rho = this->rho().field();

    //    scalarField rDeltaT(rho.size(),1/mesh_.time().deltaTValue());
    const scalarField &rDeltaT=getRDeltaT();

    PtrList<volScalarField>& Y = thermo_.Y();
    PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<scalarField> pSat(thermo_.pSat(activeSpecies));
    PtrList<scalarField> D(thermo_.diffusivity().Deff());

    const sectionalDistribution& dist = system_->distribution();
    sectionalInterpolation& interp = system_->interpolation();

    PtrList<section>& sections = system_->distribution().sections();

    const scalarField dcm(this->meanDiameter(1,0));
    const scalarField rhol(thermo_.thermoDisp().rho());

    // Nucleation

    if (nucleation_->modelType() != "none")
    {
        PtrList<scalarField> rhoDisp(thermo_.rhoDisp(activeSpecies));
        PtrList<scalarField> sigma(thermo_.sigma(activeSpecies));

        forAll(rho, celli)
        {
            const nucData ndata
            (
                nucleation_->rate
                (
                    p[celli],
                    T[celli],
                    entryList(Y,celli),
                    entryList(pSat,celli),
                    entryList(D,celli),
                    entryList(rhoDisp,celli),
                    entryList(sigma,celli)
                )
            );

            if (ndata.active())
            {
                J_.field()[celli] = ndata.J();

                interp.addToM
                (
                    ndata.s(),
                    ndata.J()/rho[celli]/rDeltaT[celli],
                    celli
                );

                const scalar Inuc(ndata.s()*ndata.J()/rho[celli]);

                forAll(activeSpecies, j)
                {
                    const scalar dZj
                    (
                        min(Inuc*ndata.z()[j], Y[j][celli])
                      / rDeltaT[celli]
                    );

                    Z[j][celli] += dZj;
                    Y[j][celli] -= dZj;
                }
            }
        }
    }

    // Condensation

    if (condensation_->modelType() != "none")
    {
        PtrList<scalarField> rhoCont(thermo_.rhoCont(contSpecies));

        forAll(rho, celli)
        {
            const conData cdata
            (
                condensation_->rate
                (
                    p[celli],
                    T[celli],
                    entryList(Y,celli),
                    entryList(Z,celli),
                    entryList(pSat,celli),
                    entryList(D,celli),
                    entryList(rhoCont,celli)
                )
            );

            if (cdata.active())
            {
                scalarList M0(dist.size(), 0.0);

                scalar sumM(0.0);

                forAll(sections, i)
                {
                    M0[i] = max(sections[i].M().field()[celli],0.0);

                    sumM += M0[i];

                    sections[i].M().field()[celli] = 0.0;
                }

                const scalarList Y0(entryList(Y,celli));
                const scalarList Z0(entryList(Z,celli));

                const scalar d(min(dcm[celli],dMax_));

                scalar dAlpha(0.0);

                forAll(activeSpecies, j)
                {
                    const scalar a((Y0[j]+Z0[j])*cdata.source()[j]);

                    const scalar b
                    (
                        max
                        (
                            (cdata.source()[j]+cdata.sink()[j]),
                            VSMALL
                        )
                    );

                    Z[j][celli] = max
                    (
                        a/b
                      + (Z0[j] - a/b)
                        * Foam::exp(-b*sumM*d/rDeltaT[celli]),
                        0.0
                    );

                    Y[j][celli] = max(Y0[j]+Z0[j]-Z[j][celli],0.0);

                    dAlpha += (Z[j][celli]-Z0[j]);

                    I_[j].field()[celli] =
                        rho[celli]*(Z[j][celli]-Z0[j])*rDeltaT[celli];
                }

                // I(s) ~ const

                const scalar Gamma
                (
                    dAlpha*rDeltaT[celli]
                  / (max(d*sumM,VSMALL))
                );

/*
                // I(s)/d ~ const

                const scalar Gamma
                (
                    dAlpha
                  * Foam::pow(6.0/(rhol[celli]*pi), 1.0/3.0)*rDeltaT[celli]
                  / (max(d*sumM,VSMALL))
                );
*/

                forAll(sections, i)
                {
                    if (M0[i] > SMALL)
                    {
                        // I(s) ~ const

                        // Note: s can become negative

                        const scalar s
                        (
                            dist[i].x()
                            + Gamma/rDeltaT[celli]*dist[i].d(rhol[celli])
                        );

/*
                        // I(s)/d ~ const

                        // Note: without clipping the argument of the 3/2 power
                        // can become negative. The clipping causes an
                        // inconsistency in the mass transfer in f

                        const scalar s
                        (
                            Foam::pow
                            (
                                max
                                (
                                    Foam::pow(dist[i].x(), 2.0/3.0)
                                    + 2.0/3.0*Gamma/rDeltaT[celli],
                                    0.0
                                ),
                                3.0/2.0
                            )
                        );
*/

                        if (s >= dist.xMin())
                        {
                            interp.addToM(s, M0[i], celli);
                        }
                    }
                }
            }
        }
    }

    // Coalescence

    if (coalescence_->modelType() != "none")
    {
        const scalarField mug(thermo_.thermoCont().mu());
        const scalarField rhog(thermo_.thermoCont().rho());

        const PtrList<coalescencePair>& pairs =
            system_->coalescencePairs();

        if (pairs.size() == 0)
        {
            system_->generateCoalescencePairs();
        }

        forAll(rho, celli)
        {
            const coaData cdata
            (
                coalescence_->rate
                (
                    p[celli],
                    T[celli],
                    mug[celli],
                    rhog[celli],
                    rhol[celli],
                    dcm[celli]
                )
            );

            if (cdata.active())
            {
                scalarList d(sections.size(), 0.0);
                scalarList M0(dist.size(), 0.0);

                forAll(sections, i)
                {
                    d[i] = dist[i].d(rhol[celli]);
                    M0[i] = max(sections[i].M().field()[celli], 0.0);
                }

                forAll(pairs, k)
                {
                    const coalescencePair& pair = pairs[k];

                    const label i(pair.i());
                    const label j(pair.j());

                    scalar beta(0.0);

                    forAll(cdata.w(), l)
                    {
                        beta +=
                            cdata.w()[l]
                          * (
                                pow(d[i],cdata.p()[l])*pow(d[j],cdata.q()[l])
                              + pow(d[i],cdata.q()[l])*pow(d[j],cdata.p()[l])
                            );
                    }

                    scalar& Mi = sections[i].M().field()[celli];
                    scalar& Mj = sections[j].M().field()[celli];

                    const scalar f
                    (
                        min
                        (
                            M0[i]*M0[j]*rho[celli]*beta/rDeltaT[celli],
                            min(Mi, Mj)
                        )
                    );

                    Mi -= f;
                    Mj -= f;

                    interp.addToM(pair.idata(), pair.s(), f, celli);
                }
            }
        }
    }

    system_->rescale();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::fixedSectional::fixedSectional
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    aerosolModel(modelType, mesh, aerosolProperties),
    system_(),
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("J", dimless/dimVolume/dimTime, 0)
    ),
    I_(thermo_.activeSpecies().size())
{
    system_.set(
        new fixedSectionalSystem(*this, coeffs())
    );

    const speciesTable& activeSpecies = thermo_.activeSpecies();

    forAll(activeSpecies, j)
    {
        I_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word(activeSpecies[j]+":I"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("I", dimMass/dimTime/dimVolume, 0)
            )
        );
    }

    phiEff_.setSize(system_->distribution().size());

    forAll(system_->distribution(), i)
    {
        const section& sec = system_->distribution()[i];

        fields_.add(sec.M());

        mesh.setFluxRequired(sec.M().name());

        phiEff_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phiEff", sec.M().name()),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("phi", dimless/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::fixedSectional::~fixedSectional()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::fixedSectional::correctModel()
{
    forAll(system_->distribution(), i)
    {
        section& sec = system_->distribution()[i];

        if(!sec.validM())
        {
            sec.initM(word(coeffs().lookup("initFromPatch")));
        }
    }

    updateDrift();
}


void Foam::aerosolModels::fixedSectional::solvePre()
{}


void Foam::aerosolModels::fixedSectional::solvePost()
{
    solveSpatial();
    solveInternal();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::aerosolModels::fixedSectional::R(const volScalarField& Y) const
{
    tmp<volScalarField> I
    (
        new volScalarField
        (
            IOobject
            (
                Y.name() + ":I:zero",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("I", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    return fvm::Su(I, Y);
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::Qdot() const
{
    return condensation_->Qdot(I_);
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::meanDiameter
(
    const scalar p,
    const scalar q
) const
{
    dimensionedScalar dMin("d", dimLength, dMin_);
    dimensionedScalar dMax("d", dimLength, dMax_);

    return max(min(system_->meanDiameter(p,q),dMax),dMin);
}

Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::medianDiameter
(
    const scalar p
) const
{
    dimensionedScalar dMin("d", dimLength, dMin_);
    dimensionedScalar dMax("d", dimLength, dMax_);

    return max(min(system_->medianDiameter(p),dMax),dMin);
}

void Foam::aerosolModels::fixedSectional::clearRates()
{
    J_ *= 0.0;

    forAll(I_, j)
    {
        I_[j] *= 0.0;
    }
}

bool Foam::aerosolModels::fixedSectional::read()
{
    if (aerosolModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
