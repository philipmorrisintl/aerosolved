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

\file saturatedMixtureFvPatchScalarField.H
\brief Saturated mixture BC

This BC can be used to specify a perfectly saturated mixture at a patch.
Alternatively, by specifying saturation ratios for each species, the mixture can
be made non-saturated. The composition of the inert gas should be configured
using the 'inertMoleFrac' entry.

*/

#ifndef saturatedMixtureFvPatchScalarField_H
#define saturatedMixtureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "aerosolModel.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class saturatedMixtureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class saturatedMixtureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data

        //- Switch if fields are checked
        Switch checkFields_;

        //- List of saturation fractions (default to 1)
        HashTable<scalar> S_;

        //- List of inert mole fractions
        HashTable<scalar> f_;


    // Private Member Functions

        //- Check if all fields have this BC
        void checkFields();

public:

    //- Runtime type information
    TypeName("saturatedMixture");


    // Constructors

        //- Construct from patch and internal field
        saturatedMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        saturatedMixtureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given saturatedMixtureFvPatchScalarField
        // onto a new patch
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new saturatedMixtureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        saturatedMixtureFvPatchScalarField
        (
            const saturatedMixtureFvPatchScalarField&,
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
                new saturatedMixtureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
