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

#include "sectionalList.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sectionalList, 0);
addToRunTimeSelectionTable(sectionalDistribution, sectionalList, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void sectionalList::update()
{
    x_ = dict().lookupType<scalarList>("x");
    y_ = dict().lookupType<scalarList>("y");

    if (x_.size() < 1)
    {
        FatalErrorInFunction
            << "Minimum sectional distribution size is 1"
            << abort(FatalError);
    }

    if (y_.size() != (x_.size()+1))
    {
        FatalErrorInFunction
            << "The y list should have size (size(x)+1). Current size of x is "
            << x_.size() << " and y is " << y_.size() << abort(FatalError);
    }

    setSize(x_.size());

    forAll(x_, i)
    {
        set
        (
            i,
            new section
            (
                i,
                sectionName(i+1),
                x_[i],
                y_[i],
                y_[i+1],
                this->sizeDimensions(),
                aerosol_.mesh()
            )
        );
    }
}

label sectionalList::search(const scalar& s, const scalarList& x) const
{
    int imin = 0;
    int imax = size()-1;
    int imid = (imax+imin)/2;

    while (imin != imax-1)
    {
        if (x[imid] <= s)
        {
            imin = imid;
        }
        else
        {
            imax = imid;
        }

        imid = (imax+imin)/2;
    }

    return imin;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sectionalList::sectionalList
(
    aerosolModel& aerosol,
    const dictionary& dict,
    const dimensionSet& dimensions
)
:
    sectionalDistribution(type(), aerosol, dict, dimensions)
{
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sectionalList::~sectionalList()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

label sectionalList::search(const scalar& s) const
{
    return search(s, y_);
}

label sectionalList::searchLower(const scalar& s) const
{
    return search(s, x_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
