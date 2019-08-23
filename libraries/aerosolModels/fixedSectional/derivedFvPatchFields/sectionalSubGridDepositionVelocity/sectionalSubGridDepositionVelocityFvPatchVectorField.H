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

\file sectionalSubGridDepositionVelocityFvPatchVectorField.H
\brief Sectional sub-grid deposition velocity BC

Uses the subGridDepositionModel to determine the deposition velocity at the
patch.

*/

#ifndef sectionalSubGridDepositionVelocityFvPatchVectorField_H
#define sectionalSubGridDepositionVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class sectionalSubGridDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sectionalSubGridDepositionVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private data

        //- Gravitational accelleration vector
        vector g_;

        //- Maximum number of Newton iterations
        label maxIter_;

        //- Newton iteration tolerance
        scalar tolerance_;


    // Private member functions

        //- Find the section index to which this BC's field belongs to
        label sectionIndex() const;


public:

    //- Runtime type information
    TypeName("sectionalSubGridDepositionVelocity");


    // Constructors

        //- Construct from patch and internal field
        sectionalSubGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sectionalSubGridDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch  onto a new patch
        sectionalSubGridDepositionVelocityFvPatchVectorField
        (
            const sectionalSubGridDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new sectionalSubGridDepositionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sectionalSubGridDepositionVelocityFvPatchVectorField
        (
            const sectionalSubGridDepositionVelocityFvPatchVectorField&,
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
                new sectionalSubGridDepositionVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
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
