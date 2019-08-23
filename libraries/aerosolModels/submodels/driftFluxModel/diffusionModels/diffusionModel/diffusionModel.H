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

\file diffusionModel.H
\brief Base class of the diffusion drift model

*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "PtrList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        diffusionModel(const diffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionModel&);


public:

    //- Runtime type information
    TypeName("diffusionModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             diffusionModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        diffusionModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected diffusion model
        static autoPtr<diffusionModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~diffusionModel();


    // Member Functions

        //- Compute the diffusivity given a specie index
        virtual tmp<volScalarField> D(const label j) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
