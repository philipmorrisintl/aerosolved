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
#include "twoMomentLogNormalAnalytical.H"
#include "fv.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    defineTypeNameAndDebug(twoMomentLogNormalAnalytical, 0);
    addToRunTimeSelectionTable
    (
        aerosolModel,
        twoMomentLogNormalAnalytical,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormalAnalytical::twoMomentLogNormalAnalytical
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    twoMomentLogNormal(modelType, mesh, aerosolProperties)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::twoMomentLogNormalAnalytical::
~twoMomentLogNormalAnalytical()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::twoMomentLogNormalAnalytical::correctModel()
{
    twoMomentLogNormal::updateDrift();
}

void Foam::aerosolModels::twoMomentLogNormalAnalytical::solvePre()
{
    twoMomentLogNormal::solvePre();
}

void Foam::aerosolModels::twoMomentLogNormalAnalytical::solvePost()
{
    // TODO: check if this call is necessary or redundant with a previous call
    clearRates();

    twoMomentLogNormal::solvePost();

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

    const scalarField &rDeltaT=getRDeltaT();
    //    scalarField rDeltaT(rho.size(),1/mesh_.time().deltaTValue());

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

                M[celli] += ndata.J()/rho[celli]/rDeltaT[celli];

                const scalar Inuc(ndata.s()*ndata.J()/rho[celli]);

                forAll(activeSpecies, j)
                {
                    Z[j][celli] +=
                        min(Inuc*ndata.z()[j], Y[j][celli])/rDeltaT[celli];

                    Y[j][celli] -=
                        min(Inuc*ndata.z()[j], Y[j][celli])/rDeltaT[celli];
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
                    const scalar Y0(Y[j][celli]);
                    const scalar Z0(Z[j][celli]);

                    const scalar Icoeffj
                    (
                        cdata.source()[j]*Y0
                      - cdata.sink()[j]*Z0
                    );

                    Icoeff += Icoeffj;

                    const scalar a((Y0+Z0)*cdata.source()[j]);

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
                      + (Z0 - a/b)
                      * Foam::exp(-b*M[celli]*dcm[celli]/rDeltaT[celli]),
                        0.0
                    );

                    Y[j][celli] = max(Y0+Z0-Z[j][celli],0.0);

                    I_[j].field()[celli] =
                        rho[celli]*(Z[j][celli]-Z0)*rDeltaT[celli];
                }

                // Complete evaporation

                if (Icoeff < 0.0)
                {
                    const scalar M0(M[celli]);

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

                    M[celli] = M0*Foam::exp(f*Icoeff*dMin_/rDeltaT[celli]);

                    J_.field()[celli] +=
                        rho[celli]*(M[celli]-M0)*rDeltaT[celli];
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
                scalar f(0.0);

                forAll(cdata.w(), i)
                {
                    f +=
                        cdata.w()[i]
                      * pow(CMD[celli], cdata.p()[i]+cdata.q()[i])
                      * exp
                        (
                            0.5
                          * (sqr(cdata.p()[i])+sqr(cdata.q()[i]))
                          * sigmaSqr
                        );
                }

                M[celli] /= (1+rho[celli]*M[celli]*f/rDeltaT[celli]);
            }
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::aerosolModels::twoMomentLogNormalAnalytical::R
(
    const volScalarField& Y
) const
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
Foam::aerosolModels::twoMomentLogNormalAnalytical::Qdot() const
{
    return twoMomentLogNormal::Qdot();
}


bool Foam::aerosolModels::twoMomentLogNormalAnalytical::read()
{
    if (twoMomentLogNormal::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
