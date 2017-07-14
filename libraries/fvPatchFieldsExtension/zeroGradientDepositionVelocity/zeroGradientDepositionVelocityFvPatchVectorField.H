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

\file zeroGradientDepositionVelocityFvPatchVectorField.H
\brief Simple zero-gradient BC for a positive deposition velocity

This BC is described in Sec. 7.2.4 in Eq. (7.23) and (7.24) in \cite thesis. If
the cell-centered value of the velocity has a wall-ward component, then the face
value is set equal to this inward component, allowing for non-zero droplet flux
accross the face and driving deposition. If the cell-centered value of the
velocity has a negative wall-ward component, then the face value of the velocity
is set to zero. This means that as long as the velocity is 'positive', i.e.,
directed at the wall, then this BC is essentially a zero-gradient BC.

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
    // Private data


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

        //- Construct by mapping given zeroGradientDepositionVelocityFvPatchVectorField
        //- onto a new patch
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

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
