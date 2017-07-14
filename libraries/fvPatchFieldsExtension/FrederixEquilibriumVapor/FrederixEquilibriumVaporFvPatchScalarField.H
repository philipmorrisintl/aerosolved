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

\file FrederixEquilibriumVaporFvPatchScalarField.H
\brief An inlet BC for a stable liquid-vapor mixture

This boundary condition sets all vapor species to be in equilibrium with the
droplet phase. This can be used as a stable inlet condition such that at the
inlet, at the specified temperature, no immediate nucleation, condensation or
evaporation takes place. This is as if the flow and aerosol dynamics are 'well
developed' as they enter the domain. If species are specified to be
non-phase-changing in the fluidThermo dictionary, then they must have an
'inertFraction' property specified, i.e., the fraction with which they
contribute to the inert mixture. This allows to specify air as a mixture of
inert gasses. There should be at least one phase-changing species and one
non-phase-changing species specified. Moreover, all \f$Y_j\f$-fields should
have this BC at a patch if one \f$Y_j\f$-field has this BC. Finally, the sum of
all inertFraction properties for all non-phase-changing species should sum up to
unity. Note that this BC is based on the standard inlet/outlet BC of OpenFOAM,
meaning that if the patch at which this BC is active is an outflow patch, then a
zero-gradient BC is applied.

*/

#ifndef FrederixEquilibriumVaporFvPatchScalarField_H
#define FrederixEquilibriumVaporFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class FrederixEquilibriumVaporFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class FrederixEquilibriumVaporFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Is this species non-phase changing?
        Switch phaseChange_;

        //- If non-phase changing this is the fraction with which this species
        //  contributes to the inert mixture
        scalar inertFraction_;

        //- Switch if checkFields() was succesfully executed
        Switch checkFields_;

        //- Name of the flux transporting the field
        word phiName_;

public:

    //- Runtime type information
    TypeName("FrederixEquilibriumVapor");


    // Constructors

        //- Construct from patch and internal field
        FrederixEquilibriumVaporFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        FrederixEquilibriumVaporFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given FrederixEquilibriumVaporFvPatchScalarField
        //- onto a new patch
        FrederixEquilibriumVaporFvPatchScalarField
        (
            const FrederixEquilibriumVaporFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        FrederixEquilibriumVaporFvPatchScalarField
        (
            const FrederixEquilibriumVaporFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new FrederixEquilibriumVaporFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        FrederixEquilibriumVaporFvPatchScalarField
        (
            const FrederixEquilibriumVaporFvPatchScalarField&,
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
                new FrederixEquilibriumVaporFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Function to check if all Y fields carry this boundary condition
            //- on this patch.
            virtual void checkFields();

            //- Function to find the species index inside the species dictionary
            virtual label speciesIndex();


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
