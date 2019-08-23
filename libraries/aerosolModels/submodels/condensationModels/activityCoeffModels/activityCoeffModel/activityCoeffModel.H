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

/**

\file activityCoeffModel.H
\brief Base class for the activity coefficient model

*/

#ifndef activityCoeffModel_H
#define activityCoeffModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class activityCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class activityCoeffModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        activityCoeffModel(const activityCoeffModel&);

        //- Disallow default bitwise assignment
        void operator=(const activityCoeffModel&);


public:

    //- Runtime type information
    TypeName("activityCoeffModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             activityCoeffModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        activityCoeffModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected activity coeff model
        static autoPtr<activityCoeffModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~activityCoeffModel();


    // Member Functions

        // Evolution

            //- Compute the activity coeff data
            virtual scalarList activity(const scalarList& Z) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
