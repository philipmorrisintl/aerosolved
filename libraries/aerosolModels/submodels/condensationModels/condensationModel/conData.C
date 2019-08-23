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

#include "conData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

conData::conData(const conData& data)
:
    source_(data.source_),
    sink_(data.sink_),
    active_(data.active_)
{}

conData::conData(const label N)
:
    source_(N, 0.0),
    sink_(N, 0.0),
    active_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

conData::~conData()
{}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void conData::operator*=(const scalar& s)
{
    forAll(source_, i)
    {
        source_[i] *= s;
        sink_[i] *= s;
    }
}

void conData::operator/=(const scalar& s)
{
    forAll(source_, i)
    {
        source_[i] /= s;
        sink_[i] /= s;
    }
}

void conData::operator=(const conData& d)
{
    source_ = d.source();
    sink_ = d.sink();
    active_ = d.active();
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

conData operator*
(
    const conData& d,
    const scalar& s
)
{
    conData data(d);

    forAll(d.source(), i)
    {
        data.source()[i] *= s;
        data.sink()[i] *= s;
    }

    return data;
}

conData operator*
(
    const scalar& s,
    const conData& d
)
{
    return d*s;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
