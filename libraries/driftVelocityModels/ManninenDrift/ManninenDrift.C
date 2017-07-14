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

#include "driftVelocityModel.H"
#include "ManninenDrift.H"
#include "fvc.H"

#include "makeDriftVelocityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{
    makeDriftVelocityTypes(ManninenDrift, driftVelocityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModels::ManninenDrift::ManninenDrift
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    driftVelocityModel(mesh, aerosol),
    g_("g", dimVelocity/dimTime, vector(0, 0, 0))
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModels::ManninenDrift::~ManninenDrift()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::driftVelocityModels::ManninenDrift::updateDropletVelocity
(
    Foam::volVectorField& V,
    const Foam::volScalarField& d,
    const Foam::volScalarField& D,
    const Foam::volVectorField& G
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    surfaceScalarField phi(fvc::interpolate(U) & mesh_.Sf());

    volVectorField dUdt
    (
        fvc::ddt(U)
      + fvc::div(phi, U, "div(phi,V)")
      - fvc::Sp(fvc::div(phi), U)
    );

    V = U - (dUdt - G) / D;

    V.correctBoundaryConditions();
}

void Foam::driftVelocityModels::ManninenDrift::updateDropDriftVelFields()
{
    // Must be a sectional aerosol model

    if (aerosol_.modType() != SECTIONALAEROSOLMODEL)
    {
        FatalErrorIn("Foam::driftVelocityModels::ManninenDrift::getDriftVelFields()")
            << "Must be a sectional model" << exit(FatalError);
    }

    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar k("k", dimEnergy/dimTemperature, 1.3806488E-23);

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");
    const volScalarField& T = thermo_.T();

    tmp<volScalarField> tRhol = thermo_.rhoLiquid();
    volScalarField& rhol = tRhol();

    tmp<volScalarField> tRhov = thermo_.rhoVapor();
    volScalarField& rhov = tRhov();

    thermo().limitLiquidDensity(rhol);
    thermo().limitVaporDensity(rhov);

    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const tmp<volScalarField> mg = thermo_.mVapor();

    // Mean free path

    const volScalarField lambda
    (
        sqrt(8.0 * k * T / (pi * mg)) * (4.0/5.0 * muEff / (p0 + p1))
    );

    forAll(aerosol_.x(), i)
    {
        const dimensionedScalar m("z", dimMass, aerosol_.x()[i]);

        const volScalarField d(pow(m/rhol * 6.0/pi, 1.0/3.0));

        const volScalarField Kn(lambda/d);

        volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

        if (!Cunningham_)
        {
            C == dimensionedScalar("one", dimless, 1.0);
        }

        const volScalarField D(18.0 * muEff / (sqr(d) * rhol * C));

        const volVectorField G((rhol-rhov)/rhol * g_);

        updateDropletVelocity(aerosol_.V()[i], d, D, G);
    }
}

void Foam::driftVelocityModels::ManninenDrift::updateDropDriftVelField
(
    const Foam::volScalarField& d
)
{
    // Must be a moment aerosol model

    if (aerosol_.modType() != MOMENTAEROSOLMODEL)
    {
        FatalErrorIn("Foam::driftVelocityModels::ManninenDrift::getDriftVelFields()")
            << "Must be a moment model" << exit(FatalError);
    }

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");
    const volScalarField& T = thermo_.T();

    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar k("k", dimEnergy/dimTemperature, 1.3806488E-23);

    tmp<volScalarField> tRhol = thermo_.rhoLiquid();
    volScalarField& rhol = tRhol();

    tmp<volScalarField> tRhov = thermo_.rhoVapor();
    volScalarField& rhov = tRhov();

    thermo().limitLiquidDensity(rhol);
    thermo().limitVaporDensity(rhov);

    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const tmp<volScalarField> mg = thermo_.mVapor();

    // Mean free path

    const volScalarField lambda
    (
        sqrt(8.0 * k * T / (pi * mg)) * (4.0/5.0 * muEff / (p0 + p1))
    );

    const dimensionedScalar smalld("d", dimLength, 1E-10);

    const volScalarField Kn(lambda/d);

    volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

    if (!Cunningham_)
    {
        C == dimensionedScalar("one", dimless, 1.0);
    }

    const volScalarField D(18.0 * muEff / (sqr(stabilise(d, smalld)) * rhol * C));

    const volVectorField G((rhol-rhov)/rhol * g_);

    updateDropletVelocity(aerosol_.V()[0], d, D, G);
}

bool Foam::driftVelocityModels::ManninenDrift::read()
{
    if (driftVelocityModel::read())
    {
        params_.lookup("g") >> g_.value();
        params_.lookup("Cunningham") >> Cunningham_;

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
