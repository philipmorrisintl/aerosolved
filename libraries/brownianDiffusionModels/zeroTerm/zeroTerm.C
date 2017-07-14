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
#include "makeBrownianDiffusionTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace brownianDiffusionModels
{
    makeDiffusionModelTypes(zeroTerm, brownianDiffusionModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::brownianDiffusionModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    brownianDiffusionModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::brownianDiffusionModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::brownianDiffusionModels::zeroTerm::update()
{
}

bool Foam::brownianDiffusionModels::zeroTerm::read()
{
    return true;
}


// ************************************************************************* //
