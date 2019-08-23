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

\file coalescenceModel.H
\brief Base class of the coalescence model

The coalescence model provides the coalescence rate function in the form of the
coaData object. This coaData object (which is essentially a polynomial) can be
evaluated by the model to determine the coalescence rate of two particles. The
model is designed in such a way that it can interact with both sectional and
moment models.

*/

#ifndef coalescenceModel_H
#define coalescenceModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "coaData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coalescenceModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        coalescenceModel(const coalescenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceModel&);


public:

    //- Runtime type information
    TypeName("coalescenceModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             coalescenceModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        coalescenceModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected coalescence model
        static autoPtr<coalescenceModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceModel();


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
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
