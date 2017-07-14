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

\file twoMomentLogNormalFvPatchScalarField.H
\brief Provide a Dirichlet BC for M based on a log-normal assumption.

This BC can be used to set a droplet number concentration based on a specified
sigma and count median diameter (CMD), such that the first moment of the
resulting distribution is equil to the total mass concentration inside the
domain as given by the inlet liquid mass fractions \f$\mathbf{Z}\f$. Note that
this BC is based on the standard inlet/outlet BC of OpenFOAM, meaning that if
the patch at which this BC is active is an outflow patch, then a zero-gradient
BC is applied.

*/

#ifndef twoMomentLogNormalFvPatchScalarField_H
#define twoMomentLogNormalFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class twoMomentLogNormalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class twoMomentLogNormalFvPatchScalarField
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
        //- onto a new patch
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

        // Access

            //- Give full name of current fv field
            inline word name();

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
