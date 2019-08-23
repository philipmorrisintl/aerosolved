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

#include "coalescencePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coalescencePair::coalescencePair(const coalescencePair& data)
:
    i_(data.i_),
    j_(data.j_),
    s_(data.s_),
    idata_(data.idata_)
{}

coalescencePair::coalescencePair
(
    const label i,
    const label j,
    const scalar s,
    const secIntData& idata
)
:
    i_(i),
    j_(j),
    s_(s),
    idata_(idata)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coalescencePair::~coalescencePair()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
