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

#ifndef Zhang_H
#define Zhang_H

#include "activityCoeffModel.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Zhang Declaration
\*---------------------------------------------------------------------------*/

class Zhang
:
    public activityCoeffModel
{
    // Private Data

        //- Name of the first specie
        word firstSpecieName_;

        //- Constant coefficients
        scalar C1_;
        scalar C2_;


    // Private member functions

        //- Disallow default bitwise copy construct
        Zhang(const Zhang&);

        //- Disallow default bitwise assignment
        void operator=(const Zhang&);


public:

    //- Runtime type information
    TypeName("Zhang");


    // Constructors

        //- Construct from aerosol model
        Zhang(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~Zhang();


    // Member Functions

        // Evolution

            //- Compute the activity coeff data
            virtual scalarList activity(const scalarList& Z) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
