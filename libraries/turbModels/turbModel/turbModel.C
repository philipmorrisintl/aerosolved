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

#include "turbModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(turbModel, 0);
    defineRunTimeSelectionTable(turbModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbModel::turbModel
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "turbulenceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    coeffs_(subDict("Coeffs")),
    velocityName_(coeffs_.lookup(word("velocityName"))),
    densityName_(coeffs_.lookup(word("densityName"))),
    massFluxName_(coeffs_.lookup(word("massFluxName"))),
    viscosityName_(coeffs_.lookup(word("viscosityName"))),
    diffusivityName_(coeffs_.lookup(word("diffusivityName"))),
    muTurb_
    (
        IOobject
        (
            "muTurb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("muTurb", dimMass/dimLength/dimTime, 0.0)
    ),
    kTurb_
    (
        IOobject
        (
            "kTurb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kTurb", dimPower/dimLength/dimTemperature, 0.0)
    ),
    U_(mesh_.lookupObject<volVectorField>( velocityName_ )),
    rho_(mesh_.lookupObject<volScalarField>( densityName_ )),
    phi_(mesh_.lookupObject<surfaceScalarField>( massFluxName_ )),
    muEff_(mesh_.lookupObject<volScalarField>( viscosityName_ )),
    kEff_(mesh_.lookupObject<volScalarField>( diffusivityName_ )),
    delta_
    (
        IOobject
        (
            "delta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("delta", dimLength, SMALL),
        "zeroGradient"
    ),
    deltaCoeff_(1.0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbModel::~turbModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

namespace Foam
{
tmp<fvVectorMatrix> Foam::turbModel::divDevRhoBeff(volVectorField& U) const
{
    return
    (
       fvm::laplacian((muEff()+muTurb()), U) + fvc::div((muEff()+muTurb())*dev2(T(fvc::grad(U))))
    );
}

}

// ************************************************************************* //
