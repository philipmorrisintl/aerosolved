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
#include "makeDiffusionTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
    makeDiffusionModelTypes(zeroTerm, diffusionModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    diffusionModel(mesh, thermo)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionModels::zeroTerm::update()
{
    forAll(thermo_.species(), j)
    {
        DY_[j] *= 0;
    }

    phic_ *= 0;
}

bool Foam::diffusionModels::zeroTerm::read()
{
    return true;
}


// ************************************************************************* //
