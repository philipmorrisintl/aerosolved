/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2017 Philip Morris International

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

\file saturatedVaporMixtureFvPatchScalarField.H
\brief An inlet BC for a stable vapors-only mixture

This BC sets the vapor mass fraction to be at saturation level. It is based on
an inlet-outlet BC so that when the sign of the flux changes this BC reverts to
a zeroGradient BC.

*/

#ifndef saturatedVaporMixtureFvPatchScalarField_H
#define saturatedVaporMixtureFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class saturatedVaporMixtureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class saturatedVaporMixtureFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Switch if checkFields() was succesfully executed
        Switch checkFields_;

        //- Name of the flux transporting the field
        word phiName_;

public:

    //- Runtime type information
    TypeName("saturatedVaporMixture");


    // Constructors

        //- Construct from patch and internal field
        saturatedVaporMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        saturatedVaporMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given saturatedVaporMixtureFvPatchScalarField
        //- onto a new patch
        saturatedVaporMixtureFvPatchScalarField
        (
            const saturatedVaporMixtureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        saturatedVaporMixtureFvPatchScalarField
        (
            const saturatedVaporMixtureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new saturatedVaporMixtureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        saturatedVaporMixtureFvPatchScalarField
        (
            const saturatedVaporMixtureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new saturatedVaporMixtureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Function to check if all Y fields carry this boundary condition
            //- on this patch.
            virtual void checkFields();

            //- Function to find the species index inside the species dictionary
            virtual label speciesIndex();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
