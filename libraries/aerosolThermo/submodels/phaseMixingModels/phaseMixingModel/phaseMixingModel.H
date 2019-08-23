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

\file phaseMixingModel.H
\brief Base class for the phase mixing model

The phaseMixingModel object determines how a property of both the continuous
and the dispersed phase is mixed into a mixture property. The mixing can be
done for a volField, cell list or patch.

*/

#ifndef phaseMixingModel_H
#define phaseMixingModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward decleration
class aerosolThermo;

/*---------------------------------------------------------------------------*\
                       Class phaseMixingModel Declaration
\*---------------------------------------------------------------------------*/

class phaseMixingModel
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseMixingModel(const phaseMixingModel&);

        //- Disallow default bitwise assignment
        void operator=(const phaseMixingModel&);


protected:

    // Protected Data

        //- Reference to the aeroso thermo object
        const aerosolThermo& thermo_;


public:

    //- Runtime type information
    TypeName("phaseMixingModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             phaseMixingModel,
             dictionary,
             (
                const aerosolThermo& thermo
             ),
             (thermo)
         );

    // Constructors

        //- Construct from type name, dictionary and thermo model
        phaseMixingModel
        (
            const word& modelType,
            const aerosolThermo& thermo
        );


    // Selectors

        //- Return a reference to the selected mixing model
        static autoPtr<phaseMixingModel> New
        (
            const aerosolThermo& thermo,
            const word& modelType
        );


    //- Destructor
    virtual ~phaseMixingModel();


    // Member Functions

        // Access

            //- Reference to the aerosol thermo object
            inline const aerosolThermo& thermo()
            {
                return thermo_;
            }

        // Mixing

            //- Compute the mixing for
            virtual tmp<volScalarField> mix
            (
                const volScalarField& phiCont,
                const volScalarField& phiDisp
            ) const = 0;

            //- Compute the mixing for a set of cells
            virtual tmp<scalarField> mix
            (
                const scalarField& phiCont,
                const scalarField& phiDisp,
                const labelList& cells
            ) const = 0;

            //- Compute the mixing for a patch
            virtual tmp<scalarField> mix
            (
                const scalarField& phiCont,
                const scalarField& phiDisp,
                const label patchi
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
