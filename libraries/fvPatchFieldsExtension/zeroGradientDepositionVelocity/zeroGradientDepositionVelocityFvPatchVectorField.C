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

#include "zeroGradientDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF)
{}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const zeroGradientDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::zeroGradientDepositionVelocityFvPatchVectorField::
zeroGradientDepositionVelocityFvPatchVectorField
(
    const zeroGradientDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(fcvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroGradientDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    const vectorField n(patch().nf());

    vectorField V(*this);

    forAll(V, i)
    {
        if ((V[i] & n[i]) > 0.0)
        {
            V[i] = (V[i] & n[i]) * n[i];
        }
        else
        {
            V[i] = vector::zero;
        }
    }

    operator==(V);
}

void Foam::zeroGradientDepositionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        zeroGradientDepositionVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
