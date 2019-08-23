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

#ifndef sectionalList_H
#define sectionalList_H

#include "sectionalDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sectionalList Declaration
\*---------------------------------------------------------------------------*/

class sectionalList
:
    public sectionalDistribution
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        sectionalList(const sectionalList&);

        //- Disallow default bitwise assignment
        void operator=(const sectionalList&);

        //- Update the distribution
        void update();

        //- Binary search function
        virtual label search(const scalar& s, const scalarList& x) const;


public:

    //- Runtime type information
    TypeName("list");


    // Constructors

        //- Construct from aerosol model
        sectionalList
        (
            aerosolModel& aerosol,
            const dictionary& dict,
            const dimensionSet& dimensions
        );


    //- Destructor
    virtual ~sectionalList();


    //- Member Functions

        //- Search section
        virtual label search(const scalar& s) const;

        //- Search lower section
        virtual label searchLower(const scalar& s) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
