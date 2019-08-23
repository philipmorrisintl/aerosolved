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

#include "driftFluxModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

driftFluxModel::driftFluxModel
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase
    (
        aerosol,
        dict.subDict("driftFluxModel"),
        "driftFluxModel",
        "driftFluxModel"
    ),
    diffusion_(),
    Brownian_(),
    inertial_()
{
    diffusion_ =
        diffusionModel::New(aerosol, this->dict().subDict("diffusion"));

    Brownian_ =
        BrownianModel::New(aerosol, this->dict().subDict("Brownian"));

    inertial_ =
        inertialModel::New(aerosol, this->dict().subDict("inertial"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

driftFluxModel::~driftFluxModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
driftFluxModel::phiDiffusion
(
    const volScalarField& Y,
    const volScalarField& D
) const
{
    const volScalarField& rho = aerosol_.rho();

    return
    (
      - linearInterpolate(rho*D)
      * fvc::snGrad(Y)*aerosol_.mesh().magSf()
    );
}

tmp<volVectorField> driftFluxModel::VDiffusion
(
    const volScalarField& Y,
    const volScalarField& D
) const
{
    return -D*fvc::grad(Y);
}

tmp<surfaceScalarField> driftFluxModel::phi
(
    const surfaceScalarField& phiInertial,
    const surfaceScalarField& phiBrownian,
    const volScalarField& DDisp,
    const PtrList<volScalarField>& DCont
) const
{
    tmp<surfaceScalarField> tphi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiDrift",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedScalar("phiDrift", dimMass/dimTime, 0)
        )
    );

    surfaceScalarField& phi = tphi.ref();

    const speciesTable& contSpecies = aerosol_.thermo().contSpecies();
    const speciesTable& dispSpecies = aerosol_.thermo().dispSpecies();

    forAll(contSpecies, j)
    {
        const volScalarField& Y = aerosol_.thermo().Y()[j];

        phi += phiDiffusion(Y,DCont[j]);
    }

    tmp<fv::convectionScheme<scalar>> mvConvInertial
    (
        fv::convectionScheme<Foam::scalar>::New
        (
            aerosol_.mesh(),
            aerosol_.thermo().fieldsZ(),
            phiInertial,
            aerosol_.mesh().divScheme("div(mvConv)")
        )
    );

    tmp<fv::convectionScheme<scalar>> mvConvBrownian
    (
        fv::convectionScheme<Foam::scalar>::New
        (
            aerosol_.mesh(),
            aerosol_.thermo().fieldsZ(),
            phiBrownian,
            aerosol_.mesh().divScheme("div(mvConv)")
        )
    );

    forAll(dispSpecies, j)
    {
        const volScalarField& Z = aerosol_.thermo().Z()[j];

        phi +=
            mvConvInertial->flux(phiInertial,Z)
          + mvConvBrownian->flux(phiBrownian,Z)
          + phiDiffusion(Z,DDisp);
    }

    return tphi;
}

tmp<volVectorField> driftFluxModel::V
(
    const volVectorField& VInertial,
    const volVectorField& VBrownian,
    const volScalarField& DDisp,
    const PtrList<volScalarField>& DCont
) const
{
    tmp<volVectorField> tV
    (
        new volVectorField
        (
            IOobject
            (
                "VDrift",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedVector("VDrift", dimLength/dimTime, vector::zero)
        )
    );

    volVectorField& V = tV.ref();

    const speciesTable& contSpecies = aerosol_.thermo().contSpecies();
    const speciesTable& dispSpecies = aerosol_.thermo().dispSpecies();

    forAll(contSpecies, j)
    {
        const volScalarField& Y = aerosol_.thermo().Y()[j];

        V += VDiffusion(Y,DCont[j]);
    }

    forAll(dispSpecies, j)
    {
        const volScalarField& Z = aerosol_.thermo().Z()[j];

        V += (VInertial+VBrownian)*Z + VDiffusion(Z,DDisp);
    }

    return tV;
}

tmp<volSymmTensorField> driftFluxModel::tau
(
    const surfaceScalarField& phiInertial,
    const surfaceScalarField& phiBrownian,
    const volScalarField& DDisp,
    const PtrList<volScalarField>& DCont
) const
{
    tmp<volSymmTensorField> ttau
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tau",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedSymmTensor
            (
                "tau",
                dimMass/sqr(dimTime)/dimLength,
                symmTensor::zero
            )
        )
    );

    volSymmTensorField& tau = ttau.ref();

    const volScalarField& rho = aerosol_.rho();

    const volVectorField VInertial(fvc::reconstruct(phiInertial)/rho);
    const volVectorField VBrownian(fvc::reconstruct(phiBrownian)/rho);

    const volVectorField VDrift(V(VInertial, VBrownian, DDisp, DCont));

    const speciesTable& contSpecies = aerosol_.thermo().contSpecies();
    const speciesTable& dispSpecies = aerosol_.thermo().dispSpecies();

    forAll(contSpecies, j)
    {
        const volScalarField& Y = aerosol_.thermo().Y()[j];

        const volScalarField Ys
        (
            max(fvc::average(linearInterpolate(Y)),SMALL)
        );

        const volVectorField W(VDiffusion(Y,DCont[j])/Ys - VDrift);

        tau -= rho*Y*sqr(W);
    }

    forAll(dispSpecies, j)
    {
        const volScalarField& Z = aerosol_.thermo().Z()[j];

        const volScalarField Zs
        (
            max(fvc::average(linearInterpolate(Z)),SMALL)
        );

        const volVectorField W
        (
            VInertial
          + VBrownian
          + VDiffusion(Z,DDisp)/Zs
          - VDrift
        );

        tau -= rho*Z*sqr(W);
    }

    return ttau;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
