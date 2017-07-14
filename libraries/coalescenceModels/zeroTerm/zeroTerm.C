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

#include "coalescenceModel.H"
#include "zeroTerm.H"

#include "makeCoalescenceTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    makeCoalescenceTypes(zeroTerm, coalescenceModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    coalescenceModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::coalescenceModels::zeroTerm::getCoaRateCell
(
    const Foam::scalar vi,
    const Foam::scalar vj,
    const Foam::label jCell,
    const Foam::scalar Kn
)
{
    return 0.0;
}

bool Foam::coalescenceModels::zeroTerm::read()
{
    if (coalescenceModel::read())
    {
        // Read coalescence model coefficients

        return true;
    }
    else
    {
        return false;
    }
}


