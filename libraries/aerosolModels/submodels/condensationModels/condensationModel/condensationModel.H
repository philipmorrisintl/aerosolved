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

\file condensationModel.H
\brief Base class of the condensationModel object

The condensationModel object provides the condensation rate in the form of the
conData object. The conData object provides the condensation rate coefficient,
which can be evaluated by the aerosol model do determine the condensation rate.

*/

#ifndef condensationModel_H
#define condensationModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "conData.H"
#include "activityCoeffModel.H"
#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationModel
:
    public aerosolSubModelBase
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        condensationModel(const condensationModel&);

        //- Disallow default bitwise assignment
        void operator=(const condensationModel&);


protected:

    // Protected Data

        //- Activity coefficient model
        autoPtr<activityCoeffModel> activity_;


public:

    //- Runtime type information
    TypeName("condensationModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             condensationModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        condensationModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected condensation model
        static autoPtr<condensationModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~condensationModel();


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
            ) const = 0;

            //- Heat of vaporization helper function
            virtual tmp<volScalarField> Qdot
            (
                const PtrList<volScalarField>& I
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
