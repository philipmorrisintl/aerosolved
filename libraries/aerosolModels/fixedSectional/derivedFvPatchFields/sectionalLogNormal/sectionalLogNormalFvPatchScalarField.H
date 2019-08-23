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

\file sectionalLogNormalFvPatchScalarField.H
\brief Sets a given log-normal size distribution at a patch

*/

#ifndef sectionalLogNormalFvPatchScalarField_H
#define sectionalLogNormalFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "aerosolModel.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class sectionalLogNormalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sectionalLogNormalFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Geometric standard deviation
        scalar sigma_;

        //- Moment (in length) to preserve
        scalar gamma_;

        //- The count median diameter
        autoPtr<Function1<scalar>> CMD_;


    // Private Member Functions

        //- Index of this section
        label sectionIndex() const;

        //- Log normal distribution integral
        scalarField logNormalIntegral
        (
            const scalar& xl,
            const scalar& xu,
            const scalarField& CMM
        ) const;


public:

    //- Runtime type information
    TypeName("sectionalLogNormal");


    // Constructors

        //- Construct from patch and internal field
        sectionalLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sectionalLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sectionalLogNormalFvPatchScalarField
        // onto a new patch
        sectionalLogNormalFvPatchScalarField
        (
            const sectionalLogNormalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sectionalLogNormalFvPatchScalarField
        (
            const sectionalLogNormalFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sectionalLogNormalFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sectionalLogNormalFvPatchScalarField
        (
            const sectionalLogNormalFvPatchScalarField&,
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
                new sectionalLogNormalFvPatchScalarField(*this, iF)
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
