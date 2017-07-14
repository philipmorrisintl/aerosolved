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

\file FrederixDepositionVelocityFvPatchVectorField.H
\brief The Lagrangian sub-grid model (LSM) BC for droplet velocity

This BC is based on the Lagrangian sub-grid model developed by Longest & Oldham
\cite Longest:2008aa. It uses a simplified Lagrangian equation of motion for the
droplet to find a reduced droplet velocity at the wall, assuming a linear
profile of the fluid velocity from cell centered value to zero at the wall. This
BC is discussed in Sec. 7.2.4 in \cite thesis and more elaborate in
\cite Frederix:2016ag.

Since a hybrid numerical root-finding algorithm is used a tolerance can be
specified using the TOL keyword and a maximum number of iterations by Niter. The
Cunningham correction factor can be taken into account with the Cunningham
(Switch) keyword and the gravity factor can be set using g. Note that the
inclusion of gravity makes this BC more elaborate than what is discussed in
\cite Longest:2008aa.

*/

#ifndef FrederixDepositionVelocityFvPatchVectorField_H
#define FrederixDepositionVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "Switch.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class FrederixDepositionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class FrederixDepositionVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private data

        //- Newton solver tolerance
        scalar TOL_;

        //- Newton solver max iterations
        label Niter_;

        //- Switch to enable the Cunningham correction factor
        Switch Cunningham_;

        //- Gravity vector
        vector g_;

        //- Parameter to prevent poorly-conditioned systems
        scalar SMALL_;

public:

    //- Runtime type information
    TypeName("FrederixDepositionVelocity");


    // Constructors

        //- Construct from patch and internal field
        FrederixDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        FrederixDepositionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given FrederixDepositionVelocityFvPatchVectorField
        //- onto a new patch
        FrederixDepositionVelocityFvPatchVectorField
        (
            const FrederixDepositionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new FrederixDepositionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        FrederixDepositionVelocityFvPatchVectorField
        (
            const FrederixDepositionVelocityFvPatchVectorField&,
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
                new FrederixDepositionVelocityFvPatchVectorField(*this, iF)
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

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Exponent of complex number
            inline complex expCmplx(const complex& z) const;

            //- Return true if the eigenvalues are such that there is no flow

            inline bool noFlow(const complex& l1, const complex& l2) const;

            //- Droplet trajectory solution
            inline scalar x
            (
                const scalar& t,
                const complex& l1,
                const complex& l2,
                const scalar& v,
                const scalar& g,
                const scalar& u
            ) const;

            //- Droplet velocity solution
            inline scalar dxdt
            (
                const scalar& t,
                const complex& l1,
                const complex& l2,
                const scalar& v,
                const scalar& g,
                const scalar& u
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
