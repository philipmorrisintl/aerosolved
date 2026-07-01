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

#include "Manninen.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "gaussConvectionScheme.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Manninen, 0);
addToRunTimeSelectionTable(dispersedInertialDriftModel, Manninen, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Manninen::Manninen
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    dispersedInertialDriftModel(type(), aerosol, dict)
{
    this->readBaseField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Manninen::~Manninen()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volVectorField& Manninen::V
(
    const volScalarField& d,
    const word sizeName
)
{
    volVectorField& V = this->getCachedField(sizeName);

    const fvMesh& mesh = aerosol_.mesh();

    const volScalarField& rho = aerosol_.rho();
    const volScalarField& rhol = aerosol_.thermo().thermoDisp().rho();
    const volScalarField& rhog = aerosol_.thermo().thermoCont().rho();

    const volVectorField& U = aerosol_.U();
    const surfaceScalarField& phi = aerosol_.phi();

    const volScalarField& mug = aerosol_.thermo().thermoCont().mu();
        const rhoAerosolPhaseThermo& thermoCont = aerosol_.thermo().thermoCont();
    const basicSpecieMixture& compCont = thermoCont.composition();
    const label j = thermoCont.species()[aerosol_.thermo().inertSpecie()];

    const volScalarField& p = aerosol_.thermo().p();
    const volScalarField& T = thermoCont.T();

    const scalar pi = constant::mathematical::pi;
    const scalar k = constant::physicoChemical::k.value();
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar W = compCont.W(j);
    const scalar mg = 0.001 * W / NA;

// Calculate mean free path

    volScalarField lambda
    (
        IOobject("lambda", mesh.time().timeName(), mesh),
        Foam::sqrt(pi * k * T / (2.0 * mg)) * (mug / p)
    );
    lambda.dimensions().reset(dimLength);

    const scalar dMinValue = aerosol_.dMin();
    const dimensionedScalar dMin("dMin", dimLength, dMinValue);

// Calculate Knudsen number

    volScalarField Kn
    (
        IOobject("Kn", mesh.time().timeName(), mesh),
        (2.0 * lambda) / (max(d, dMin) + dimensionedScalar("small", dimLength, SMALL))
    );
    Kn.dimensions().reset(dimless);

// Calculate Cunningham slip correction factor with safeguards

    volScalarField C
    (
        IOobject("C", mesh.time().timeName(), mesh),
        1.0 + (Kn / 2.0) * (2.34 + 1.05 * exp(-0.39 / ((Kn / 2.0) + dimensionedScalar("small", dimless, SMALL))))
    );
    C.dimensions().reset(dimless);


// Calculate relaxation time with additional safeguards


    volScalarField tau
    (
        IOobject("tau", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("tau", dimTime, 0.0)
    );

   

// Shape factor definition


        const scalar monoRad = aerosol_.monoRad();
        const scalar pfFm    = aerosol_.dysfPfFm();
        const scalar expFm   = aerosol_.dysfExpFm();
        const scalar pfTr    = aerosol_.dysfPfTr();
        const scalar expTr   = aerosol_.dysfExpTr();
        const scalar pfCont  = aerosol_.dysfPfCont();
        const scalar expCont = aerosol_.dysfExpCont();

        forAll(d, celli)
        {
            const scalar dc = max(d[celli], SMALL);
            const scalar mu = max(mug[celli], SMALL);
            const scalar c  = max(C[celli], SMALL);
            const scalar KnVal = Kn[celli];

            scalar sf = 1.0;

            if (KnVal < 0.1)
                sf = pfCont * pow(0.5 * dc / monoRad, expCont);
            else if (KnVal > 10.0)
                sf = pfFm * pow(0.5 * dc / monoRad, expFm);
            else
                sf = pfTr * pow(0.5 * dc / monoRad, expTr);

            sf = max(sf, SMALL);
            tau[celli] = sqr(dc) * rhol[celli] * c / (18.0 * mu * sf);
        }


    tau = max(tau, dimensionedScalar("minTau", dimTime, SMALL));
    
    const volVectorField G((rhol-rhog)/rhol*g_);

    const surfaceScalarField phiU(phi/linearInterpolate(rho));

    const volVectorField a
    (
        fvc::ddt(U)
      + fv::gaussConvectionScheme<vector>
        (
            mesh,
            phiU,
            upwind<vector>(mesh, phiU)
        ).fvcDiv(phiU, U)
      - fvc::Sp(fvc::div(phiU), U)
    );

    V = (G-a)*tau;
    V.correctBoundaryConditions();

    limit(V);

    return V;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
