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

\file inertialModel.H
\brief Base class of the inertial drift model

*/

#ifndef inertialModel_H
#define inertialModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class inertialModel Declaration
\*---------------------------------------------------------------------------*/

class inertialModel
:
    public aerosolSubModelBase
{
private:

    // Private Data

        //- Max relative velocity
        scalar VMax_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        inertialModel(const inertialModel&);

        //- Disallow default bitwise assignment
        void operator=(const inertialModel&);


protected:

    // Protected Data

        //- Gravitational accelleration
        uniformDimensionedVectorField g_;


public:

    //- Runtime type information
    TypeName("inertialModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             inertialModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        inertialModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected inertial model
        static autoPtr<inertialModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~inertialModel();


    // Member Functions

        //- Compute the particle Reynolds number
        tmp<volScalarField> Re
        (
            const volScalarField& d,
            const volVectorField& V
        ) const;

        //- Function to limit the velocity
        tmp<volVectorField> limit(volVectorField& U) const;

        //- Compute the relative velocity given a droplet size field and size
        //  name
        virtual tmp<volVectorField> V
        (
            const volScalarField& d,
            const word sizeName
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
