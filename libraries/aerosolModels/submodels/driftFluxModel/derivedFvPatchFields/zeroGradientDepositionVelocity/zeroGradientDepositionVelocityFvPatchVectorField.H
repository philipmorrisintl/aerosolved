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

\file zeroGradientDepositionVelocityFvPatchVectorField.H
\brief The ZG deposition BC

This BC takes the cell-centered velocity if it is wall-ward, or sets it to zero
otherwise.

*/

#ifndef zeroGradientDepositionVelocityFvPatchVectorField_H
#define zeroGradientDepositionVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class zeroGradientDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class zeroGradientDepositionVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("zeroGradientDepositionVelocity");


    // Constructors

        //- Construct from patch and internal field
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch  onto a new patch
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const zeroGradientDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new zeroGradientDepositionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroGradientDepositionVelocityFvPatchVectorField
        (
            const zeroGradientDepositionVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new zeroGradientDepositionVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
