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

\file nucleationModel.H
\brief Base class of the nucleation model

The nucleationModel object provides the aerosol model with the nucleation rate,
critical cluster size and critical cluster composition.

*/

#ifndef nucleationModel_H
#define nucleationModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "nucData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nucleationModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nucleationModel(const nucleationModel&);

        //- Disallow default bitwise assignment
        void operator=(const nucleationModel&);


public:

    //- Runtime type information
    TypeName("nucleationModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             nucleationModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        nucleationModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected nucleation model
        static autoPtr<nucleationModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~nucleationModel();


    // Member Functions

        // Evolution

            //- Compute the nucleation data
            virtual nucData rate
            (
                const scalar& p,
                const scalar& T,
                const scalarList& Y,
                const scalarList& pSat,
                const scalarList& D,
                const scalarList& rhoDisp,
                const scalarList& sigma
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
