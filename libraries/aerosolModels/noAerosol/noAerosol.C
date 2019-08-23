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
#include "noAerosol.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    defineTypeNameAndDebug(noAerosol, 0);
    addToRunTimeSelectionTable(aerosolModel, noAerosol, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::noAerosol::noAerosol
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    aerosolModel(modelType, mesh, aerosolProperties)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::noAerosol::~noAerosol()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::noAerosol::correctModel()
{}

void Foam::aerosolModels::noAerosol::solvePre()
{}

void Foam::aerosolModels::noAerosol::solvePost()
{}

Foam::tmp<Foam::fvScalarMatrix>
Foam::aerosolModels::noAerosol::R(const volScalarField& Y) const
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
Foam::aerosolModels::noAerosol::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Qdot:zero",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    return tQdot;
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::noAerosol::meanDiameter
(
    const scalar p,
    const scalar q
) const
{
    tmp<volScalarField> td
    (
        new volScalarField
        (
            IOobject
            (
                "d",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("D", dimLength, 0.0)
        )
    );

    return td;
}

Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::noAerosol::medianDiameter
(
    const scalar p
) const
{
    tmp<volScalarField> td
    (
        new volScalarField
        (
            IOobject
            (
                "d",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("D", dimLength, 0.0)
        )
    );

    return td;
}

bool Foam::aerosolModels::noAerosol::read()
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
