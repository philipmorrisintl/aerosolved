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

#include "coaData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coaData::coaData(const coaData& data)
:
    w_(data.w_),
    p_(data.p_),
    q_(data.q_),
    active_(data.active_)
{}

coaData::coaData(const label N)
:
    w_(N, 0.0),
    p_(N, 0.0),
    q_(N, 0.0),
    active_(false)
{}

coaData::coaData
(
    const scalarList& w,
    const scalarList& p,
    const scalarList& q,
    const Switch& active
)
:
    w_(w),
    p_(p),
    q_(q),
    active_(active)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coaData::~coaData()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
