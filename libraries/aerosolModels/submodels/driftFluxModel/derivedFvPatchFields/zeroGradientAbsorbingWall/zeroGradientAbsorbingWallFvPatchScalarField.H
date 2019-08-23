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

\file zeroGradientAbsorbingWallFvPatchScalarField.H
\brief The ZG absorbing wall BC

This BC uses the zero-gradient approach when evaluating the value at the patch,
or the zero-value when evaluating the gradient at the patch. This is useful to
allow for both convection and diffusion of the scalar at hand.

*/

#ifndef zeroGradientAbsorbingWallFvPatchScalarFields_H
#define zeroGradientAbsorbingWallFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class zeroGradientAbsorbingWallFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class zeroGradientAbsorbingWallFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("zeroGradientAbsorbingWall");


    // Constructors

        //- Construct from patch and internal field
        zeroGradientAbsorbingWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroGradientAbsorbingWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch onto a new patch
        zeroGradientAbsorbingWallFvPatchScalarField
        (
            const zeroGradientAbsorbingWallFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroGradientAbsorbingWallFvPatchScalarField
        (
            const zeroGradientAbsorbingWallFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroGradientAbsorbingWallFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroGradientAbsorbingWallFvPatchScalarField
        (
            const zeroGradientAbsorbingWallFvPatchScalarField&,
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
                new zeroGradientAbsorbingWallFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return the gradient at the boundary
            virtual tmp<scalarField> snGrad() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- zeroGradient coefficients
            virtual tmp<scalarField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- zeroGradient coefficients
            virtual tmp<scalarField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- zero fixedValue coefficients
            virtual tmp<scalarField> gradientInternalCoeffs() const;

            //- zero fixedValue coefficients
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
