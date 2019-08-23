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

\file twoMomentLogNormalFvPatchScalarField.H
\brief BC for the two-moment log-normal aerosol model to specify a log-normal
distribution

*/

#ifndef twoMomentLogNormalFvPatchScalarField_H
#define twoMomentLogNormalFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class twoMomentLogNormalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class twoMomentLogNormalFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- The count median diameter
        autoPtr<Function1<scalar>> CMD_;


public:

    //- Runtime type information
    TypeName("twoMomentLogNormal");


    // Constructors

        //- Construct from patch and internal field
        twoMomentLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        twoMomentLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given twoMomentLogNormalFvPatchScalarField
        // onto a new patch
        twoMomentLogNormalFvPatchScalarField
        (
            const twoMomentLogNormalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        twoMomentLogNormalFvPatchScalarField
        (
            const twoMomentLogNormalFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new twoMomentLogNormalFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        twoMomentLogNormalFvPatchScalarField
        (
            const twoMomentLogNormalFvPatchScalarField&,
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
                new twoMomentLogNormalFvPatchScalarField(*this, iF)
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
