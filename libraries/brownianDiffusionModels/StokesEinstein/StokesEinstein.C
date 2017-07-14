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

#include "StokesEinstein.H"
#include "makeBrownianDiffusionTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace brownianDiffusionModels
{
    makeDiffusionModelTypes(StokesEinstein, brownianDiffusionModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::brownianDiffusionModels::StokesEinstein::StokesEinstein
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    brownianDiffusionModel(mesh, aerosol),
    Dcrit_(0.0)
{
   read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::brownianDiffusionModels::StokesEinstein::~StokesEinstein()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::brownianDiffusionModels::StokesEinstein::update()
{
    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");
    const volScalarField& T = thermo_.T();

    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar k("k", dimEnergy/dimTemperature, 1.3806488E-23);

    tmp<volScalarField> tRhol = thermo_.rhoLiquid();
    volScalarField& rhol = tRhol();

    thermo_.limitLiquidDensity(rhol);

    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const tmp<volScalarField> mg = thermo_.mVapor();

    // Mean free path

    const volScalarField lambda
    (
        sqrt(8.0 * k * T / (pi * mg)) * (4.0/5.0 * muEff / (p0 + p1))
    );

    const dimensionedScalar Dcrit("Dcrit", dimLength, Dcrit_);

    if (aerosol_.modType() == MOMENTAEROSOLMODEL)
    {
        const volScalarField d(stabilise(aerosol_.dmm(), Dcrit));

        const volScalarField Kn(lambda/d);

        const volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

        DM_[0] = k * T * C / (3.0 * pi * muEff * d);
    }
    else if (aerosol_.modType() == SECTIONALAEROSOLMODEL)
    {
        forAll(aerosol_.x(), i)
        {
            const dimensionedScalar m("z", dimMass, aerosol_.x()[i]);

            const volScalarField d(stabilise(pow(m/rhol * 6.0/pi, 1.0/3.0), Dcrit));

            const volScalarField Kn(lambda/d);

            const volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

            DM_[i] = k * T * C / (3.0 * pi * muEff * d);
        }
    }
    else
    {
        FatalErrorIn("Foam::brownianDiffusionModels::StokesEinstein()")
            << "Invalid aerosol model type" << exit(FatalError);
    }
}

bool Foam::brownianDiffusionModels::StokesEinstein::read()
{
    coeffs_.lookup("Dcrit") >> Dcrit_;

    return true;
}


// ************************************************************************* //
