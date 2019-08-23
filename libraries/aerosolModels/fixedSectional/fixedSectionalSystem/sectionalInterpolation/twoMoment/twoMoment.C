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

#include "twoMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoMoment, 0);
addToRunTimeSelectionTable(sectionalInterpolation, twoMoment, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoMoment::twoMoment
(
    aerosolModel& aerosol,
    sectionalDistribution& distribution,
    const dictionary& dict
)
:
    sectionalInterpolation(type(), aerosol, distribution, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoMoment::~twoMoment()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

secIntData twoMoment::interp(const scalar& s) const
{
    const scalarList& x = distribution_.x();

    secIntData idata(2);

    labelList& i = idata.i();
    scalarList& w = idata.w();
    scalar& xi = idata.xi();

    i[0] = distribution_.findLower(s, true);
    i[1] = distribution_.findUpper(s, true);

    xi = max(min(s,distribution_.xMax()),distribution_.xMin());

    w[0] =  (xi-x[i[1]])/(x[i[0]]-x[i[1]]);
    w[1] = -(xi-x[i[0]])/(x[i[0]]-x[i[1]]);

    return idata;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
