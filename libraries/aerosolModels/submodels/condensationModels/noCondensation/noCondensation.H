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

#ifndef noCondensation_H
#define noCondensation_H

#include "condensationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noCondensation Declaration
\*---------------------------------------------------------------------------*/

class noCondensation
:
    public condensationModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noCondensation(const noCondensation&);

        //- Disallow default bitwise assignment
        void operator=(const noCondensation&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from aerosol model
        noCondensation(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~noCondensation();


    // Member Functions

        // Evolution

            //- Compute the condensation rate coefficients
            virtual conData rate
            (
                const scalar& p,
                const scalar& T,
                const scalarList& Y,
                const scalarList& Z,
                const scalarList& pSat,
                const scalarList& D,
                const scalarList& rhoCont
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
