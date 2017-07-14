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

\file sectionalConstantFvPatchScalarField.H
\brief Set a constant droplet number concentration for each section

This BC sets a droplet number concentration which is equal for all sections at
the patch at which this BC is active, such that the first moment of the droplet
number distribution is equal to the total mass concentration of liquid at this
patch as expressed by the liquid mass fractions \f$\mathbf{Z}\f$. Note that this
BC is based on the standard inlet/outlet BC of OpenFOAM, meaning that if the
patch at which this BC is active is an outflow patch, then a zero-gradient BC is
applied.

*/

#ifndef sectionalConstantFvPatchScalarField_H
#define sectionalConstantFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class sectionalConstantFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sectionalConstantFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;


public:

    //- Runtime type information
    TypeName("sectionalConstant");


    // Constructors

        //- Construct from patch and internal field
        sectionalConstantFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sectionalConstantFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sectionalConstantFvPatchScalarField
        //-  onto a new patch
        sectionalConstantFvPatchScalarField
        (
            const sectionalConstantFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sectionalConstantFvPatchScalarField
        (
            const sectionalConstantFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sectionalConstantFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sectionalConstantFvPatchScalarField
        (
            const sectionalConstantFvPatchScalarField&,
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
                new sectionalConstantFvPatchScalarField(*this, iF)
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
