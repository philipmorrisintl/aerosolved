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

#include "diffusionModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(diffusionModel, 0);
    defineRunTimeSelectionTable(diffusionModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModel::diffusionModel
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "diffusionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    coeffs_(subDict("Coeffs")),
    DY_(thermo.nSpecies()),
    phic_
    (
        IOobject
        (
            "phic",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phic", dimDensity*dimVelocity*dimArea, 0.0)
    )
{
    forAll(thermo.species(), i)
    {
        word speciesName = thermo.species().keys()[i];

        DY_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesName +  "DY",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("DY", dimLength*dimLength/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionModel::~diffusionModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
