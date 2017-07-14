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
#include "zeroTerm.H"

#include "makeDriftVelocityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{
    makeDriftVelocityTypes(zeroTerm, driftVelocityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    driftVelocityModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::driftVelocityModels::zeroTerm::updateDropDriftVelFields()
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    forAll(aerosol_.V(), i)
    {
        aerosol_.V()[i] == U;
        aerosol_.V()[i].correctBoundaryConditions();
    }
}

void Foam::driftVelocityModels::zeroTerm::updateDropDriftVelField
(
    const Foam::volScalarField& d
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    forAll(aerosol_.V(), i)
    {
        aerosol_.V()[i] == U;
        aerosol_.V()[i].correctBoundaryConditions();
    }
}

bool Foam::driftVelocityModels::zeroTerm::read()
{
    if (driftVelocityModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
