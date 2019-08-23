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

#ifndef noCoalescence_H
#define noCoalescence_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noCoalescence Declaration
\*---------------------------------------------------------------------------*/

class noCoalescence
:
    public coalescenceModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        noCoalescence(const noCoalescence&);

        //- Disallow default bitwise assignment
        void operator=(const noCoalescence&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from aerosol model
        noCoalescence(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~noCoalescence();


    // Member Functions

        // Evolution

            //- Compute the coalescence data
            virtual coaData rate
            (
                const scalar& p,
                const scalar& T,
                const scalar& mu,
                const scalar& rhog,
                const scalar& rhol,
                const scalar& d
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
