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

#include "zeroTerm.H"
#include "makeTurbulenceTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbModels
{
    makeTurbModelTypes(zeroTerm, turbModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    turbModel(mesh, thermo)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbModels::zeroTerm::update()
{
    muTurb_ *= 0.0;
    muTurb_.correctBoundaryConditions();
    kTurb_ *= 0.0;
    kTurb_.correctBoundaryConditions();
}


bool Foam::turbModels::zeroTerm::read()
{
    if (regIOobject::read())
    {
       return true;
    }
    else
    {
        return false;
    }

    return true;
}

Foam::tmp<Foam::volScalarField> Foam::turbModels::zeroTerm::K(volVectorField& U )
{
    volTensorField gradU = fvc::grad(U);
    return
    (
     0.0*muTurb()*mag(symm(gradU))
    );
}

// ************************************************************************* //
