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

\file sectionalLogNormalFvPatchScalarField.H
\brief Set a log-normal distribution in terms of sectional droplet number
concentrations

Using a specified sigma and count median diameter (CMD) a droplet number
distribution can be chosen in terms of the sectional droplet number
concentrations such that a log-normal distribution is mimicked (log-normal
w.r.t. droplet diameter). This BC respects the consistency relation (Eq.
(2.31)). Note that this BC is based on the standard inlet/outlet BC of OpenFOAM,
meaning that if the patch at which this BC is active is an outflow patch, then a
zero-gradient BC is applied.

*/

#ifndef sectionalLogNormalFvPatchScalarField_H
#define sectionalLogNormalFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class sectionalLogNormalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sectionalLogNormalFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Log-normal std
        scalar sigma_;

        //- Log-normal CMD
        scalar CMD_;

        //- Name of the flux transporting the field
        word phiName_;


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
        //- onto a new patch
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

        // Access

            //- Give full name of current fv field
            inline word name();

            //- Give section number
            inline label sectionNum();

            //- Check if this is a validly named field
            inline bool checkName(const word name);


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
