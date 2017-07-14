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

\file fixedMeanVelocityFvPatchVectorField.H
\brief Sets a mean velocity by pushing a plug flow or by sucking the scaled
cell-centered velocity profile

This BC can be used to specify a fixed mean velocity at an inlet our outlet, by
setting the `meanU' value in this BC, for example as follows:

\code{.unparsed}
outlet
{
    type    fixedMeanVelocity;
    meanU   0.1;
    value   uniform (0.0 0.0 0.0);
}
\endcode

This creates a mean velocity at patch 'outlet' of 0.1 m/s. In case the specified
meanU value is smaller than zero, a plug flow is pushed into the domain (the
negative sign is w.r.t. the outward normal of the patch faces). If the sign of
meanU is positive then the velocities of the cell values directly adjacent to
the patch faces are taken, copied to the patch faces and scaled such that their
surface-average is equal to the meanU value. This allows to let the outflow
develop in a natural way without forcing the flow to a certain profile at the
outlet. This prevents oscillatory and unstable solutions at the outlet. This BC
is very useful to mimic a pump at the outlet which sucks fluid from the domain.
This BC works well with the zeroGradient BC for pressure.

*/

#ifndef fixedMeanVelocityFvPatchVectorField_H
#define fixedMeanVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
            Class fixedMeanVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedMeanVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        autoPtr<DataEntry<scalar> > meanU_;

        //- Current time index
        label curTimeIndex_;


public:

   //- Runtime type information
   TypeName("fixedMeanVelocity");


   // Constructors

        //- Construct from patch and internal field
        fixedMeanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedMeanVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedMeanVelocityFvPatchVectorField
        //- onto a new patch
        fixedMeanVelocityFvPatchVectorField
        (
            const fixedMeanVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedMeanVelocityFvPatchVectorField
        (
            const fixedMeanVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedMeanVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedMeanVelocityFvPatchVectorField
        (
            const fixedMeanVelocityFvPatchVectorField&,
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
                new fixedMeanVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
