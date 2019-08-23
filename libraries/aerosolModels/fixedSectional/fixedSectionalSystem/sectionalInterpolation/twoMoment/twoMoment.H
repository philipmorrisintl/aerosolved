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

#ifndef twoMoment_H
#define twoMoment_H

#include "sectionalInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoMoment Declaration
\*---------------------------------------------------------------------------*/

class twoMoment
:
    public sectionalInterpolation
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        twoMoment(const twoMoment&);

        //- Disallow default bitwise assignment
        void operator=(const twoMoment&);


public:

    //- Runtime type information
    TypeName("twoMoment");


    // Constructors

        //- Construct from aerosol model
        twoMoment
        (
            aerosolModel& aerosol,
            sectionalDistribution& distribution,
            const dictionary& dict
        );


    //- Destructor
    virtual ~twoMoment();


    //- Member Functions

        //- Compute the indices and weights given a size
        virtual secIntData interp(const scalar& s) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
