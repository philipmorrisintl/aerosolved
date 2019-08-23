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
#include "twoMomentLogNormal.H"
#include "fv.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    defineTypeNameAndDebug(twoMomentLogNormal, 0);
    addToRunTimeSelectionTable(aerosolModel, twoMomentLogNormal, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::aerosolModels::twoMomentLogNormal::updateSources()
{
    clearRates();

    if
    (
        condensation_->modelType() == "none"
     && nucleation_->modelType() == "none"
     && coalescence_->modelType() == "none"
    )
    {
        return;
    }

    const scalar pi = constant::mathematical::pi;

    const speciesTable& activeSpecies = thermo_.activeSpecies();
    const speciesTable& contSpecies = thermo_.contSpecies();

    scalarField& M = M_.field();

    const scalarField& p = thermo_.p().field();
    const scalarField& T = thermo_.T().field();
    const scalarField& rho = this->rho().field();

    PtrList<volScalarField>& Y = thermo_.Y();
    PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<scalarField> pSat(thermo_.pSat(activeSpecies));
    PtrList<scalarField> D(thermo_.diffusivity().Deff());

    const scalarField CMD(this->medianDiameter(0));
    const scalarField rhol(thermo_.thermoDisp().rho());

    // Nucleation

    if (nucleation_->modelType() != "none")
    {
        PtrList<scalarField> rhoDisp(thermo_.rhoDisp(activeSpecies));
        PtrList<scalarField> sigma(thermo_.sigma(activeSpecies));

        forAll(M, celli)
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

                forAll(activeSpecies, j)
                {
                    I_[j].field()[celli] = ndata.s()*ndata.z()[j]*ndata.J();
                }
            }
        }
    }

    // Condensation

    if (condensation_->modelType() != "none")
    {
        const scalarField dcm(this->meanDiameter(1,0));

        PtrList<scalarField> rhoCont(thermo_.rhoCont(contSpecies));

        forAll(M, celli)
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
                scalar Icoeff(0.0);

                forAll(activeSpecies, j)
                {
                    const scalar Icoeffj
                    (
                        cdata.source()[j]*Y[j][celli]
                      - cdata.sink()[j]*Z[j][celli]
                    );

                    Icoeff += Icoeffj;

                    I_[j].field()[celli] +=
                        Icoeffj*rho[celli]*M[celli]*dcm[celli];
                }


                // Complete evaporation

                if (Icoeff < 0.0)
                {
                    const scalar sMin
                    (
                        1.0/6.0*pi
                      * Foam::pow(dMin_,3.0)*rhol[celli]
                    );

                    const scalar CMM
                    (
                        1.0/6.0*pi
                      * Foam::pow(CMD[celli],3.0)*rhol[celli]
                    );

                    const scalar f
                    (
                        1.0/(3.0*Foam::sqrt(2.0*pi)*sMin*Foam::log(sigma_))
                      * Foam::exp
                        (
                          - Foam::sqr(Foam::log(sMin/max(CMM,VSMALL)))
                          / (18.0*Foam::sqr(Foam::log(sigma_)))
                        )
                    );

                    J_.field()[celli] += M[celli]*f*Icoeff*dMin_;
                }
            }
        }
    }

    // Coalescence

    if (coalescence_->modelType() != "none")
    {
        const scalarField mug(thermo_.thermoCont().mu());
        const scalarField rhog(thermo_.thermoCont().rho());

        const scalar sigmaSqr(Foam::sqr(Foam::log(sigma_)));

        forAll(M, celli)
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
                    CMD[celli]
                )
            );

            if (cdata.active())
            {
                forAll(cdata.w(), i)
                {
                    f_[celli] +=
                        cdata.w()[i]
                      * pow(CMD[celli], cdata.p()[i]+cdata.q()[i])
                      * exp
                        (
                            0.5
                          * (sqr(cdata.p()[i])+sqr(cdata.q()[i]))
                          * sigmaSqr
                        );
                }
            }
        }
    }
}

void Foam::aerosolModels::twoMomentLogNormal::updateDrift()
{
    if (this->drift().inertial().type() != "none")
    {
        const volScalarField& rho = this->rho();

        phiInertial_ = fvc::flux
        (
            this->drift().inertial().V(this->meanDiameter(5.0,3.0)(), "d53")
          * rho
        );
    }
    else
    {
        phiInertial_ *= 0.0;
    }

    if (this->drift().Brownian().type() != "none")
    {
        DDisp_ = this->drift().Brownian().D(this->meanDiameter(3.0,2.0)());

        const volScalarField& rho = this->rho();

        phiBrownian_ = -fvc::snGrad(rho*DDisp_)*mesh_.magSf();
    }
    else
    {
        DDisp_ *= 0.0;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormal::twoMomentLogNormal
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    aerosolModel(modelType, mesh, aerosolProperties),
    M_
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
    ),
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
    f_
    (
        IOobject
        (
            "f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("f", dimVolume/dimTime, 0)
    ),
    I_(thermo_.activeSpecies().size()),
    sigma_(readScalar(coeffs_.lookup("sigma")))
{
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

    fields_.add(M_);

    mesh.setFluxRequired(M_.name());

    phiEff_.setSize(1);

    phiEff_.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("phiEff", M_.name()),
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormal::~twoMomentLogNormal()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::twoMomentLogNormal::correctModel()
{
    updateSources();
    updateDrift();
}

void Foam::aerosolModels::twoMomentLogNormal::solvePre()
{}

void Foam::aerosolModels::twoMomentLogNormal::solvePost()
{
    const volScalarField& rho = this->rho();

    const surfaceScalarField& phi = this->phi();

    const volScalarField D(turbulence().mut()+this->DDisp()*rho);

    fv::options& fvOptions(fv::options::New(mesh_));

    fvScalarMatrix MEqn
    (
        fvm::ddt(rho, M_)
      + mvPhi_->fvmDiv(phi, M_)
      + mvPhiInertial_->fvmDiv(phiInertial_, M_)
      + mvPhiBrownian_->fvmDiv(phiBrownian_, M_)
      - mvPhiDrift_->fvmDiv(phiDrift_, M_)
      ==
        fvm::laplacian(D, M_,"laplacian(D,M)")
      + fvm::Su(J_, M_)
      - fvm::Sp(f_*Foam::sqr(rho)*M_, M_)
      + fvOptions(rho, M_)
    );

    MEqn.relax();

    fvOptions.constrain(MEqn);

    MEqn.solve(mesh_.solver("M"));

    fvOptions.correct(M_);

    M_.max(0.0);

    phiEff_[0] = MEqn.flux();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::aerosolModels::twoMomentLogNormal::R(const volScalarField& Y) const
{
    const speciesTable& activeSpecies = thermo_.activeSpecies();

    const word phaseName(Y.group());
    const word speciesName(Y.member());

    if (activeSpecies.found(speciesName))
    {
        const label j(activeSpecies[speciesName]);

        if (phaseName == thermo_.phaseNameCont())
        {
            return -fvm::Su(I_[j], Y);
        }
        else if (phaseName == thermo_.phaseNameDisp())
        {
            return fvm::Su(I_[j], Y);
        }
        else
        {
            FatalErrorInFunction
                << "Phase name '" << phaseName << "' is invalid" << nl
                << abort(FatalError);

            return fvm::Su(I_[0], Y);
        }
    }
    else
    {
        return fvm::Su(I_[0]*0.0, Y);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::twoMomentLogNormal::Qdot() const
{
    return condensation_->Qdot(I_);
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::twoMomentLogNormal::meanDiameter
(
    const scalar p,
    const scalar q
) const
{
    const scalar pi = constant::mathematical::pi;

    const volScalarField alpha(thermo_.sumZ());
    const volScalarField rhol(thermo_.thermoDisp().rho());

    const dimensionedScalar smallM("M", M_.dimensions(), SMALL);

    const dimensionedScalar dMin("d", dimLength, dMin_);
    const dimensionedScalar dMax("d", dimLength, dMax_);

    return min
    (
        max
        (
            Foam::pow(alpha*6.0/(pi*rhol*max(M_,smallM)), 1.0/3.0)
          * Foam::exp(0.5*(p+q-3.0)*Foam::sqr(Foam::log(sigma_))),
            dMin
        ),
        dMax
    );
}

Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::twoMomentLogNormal::medianDiameter
(
    const scalar p
) const
{
    return meanDiameter(2.0*p, 0.0);
}

void Foam::aerosolModels::twoMomentLogNormal::clearRates()
{
    J_ *= 0.0;
    f_ *= 0.0;

    forAll(I_, j)
    {
        I_[j] *= 0.0;
    }
}

bool Foam::aerosolModels::twoMomentLogNormal::read()
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
