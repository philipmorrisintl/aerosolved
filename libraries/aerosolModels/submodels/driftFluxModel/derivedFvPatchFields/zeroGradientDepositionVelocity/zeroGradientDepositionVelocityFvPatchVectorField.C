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
    DebugInFunction << endl;

    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    const vectorField n(patch().nf());
    const vectorField c(patch().Cf());

    vectorField V(*this);

    forAll(V, facei)
    {
        if ((V[facei] & n[facei]) < 0.0)
        {
            V[facei] = V[facei] - (V[facei]&n[facei])*n[facei];
        }
    }

    operator==(V);
}


void Foam::zeroGradientDepositionVelocityFvPatchVectorField::write
(
    Ostream& os
) const
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
