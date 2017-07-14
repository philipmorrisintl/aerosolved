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

#include "conductivityModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(conductivityModel, 0);
    defineRunTimeSelectionTable(conductivityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityModel::conductivityModel
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            "conductivityProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    coeffs_(subDict("Coeffs")),
    kEff_
    (
        IOobject
        (
            "kEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kEff", dimPower/dimLength/dimTemperature, 0.0)
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityModel::~conductivityModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
