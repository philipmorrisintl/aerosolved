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

#include "zeroGradientAbsorbingWallFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroGradientAbsorbingWallFvPatchScalarField::
zeroGradientAbsorbingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


Foam::zeroGradientAbsorbingWallFvPatchScalarField::
zeroGradientAbsorbingWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict)
{}


Foam::zeroGradientAbsorbingWallFvPatchScalarField::
zeroGradientAbsorbingWallFvPatchScalarField
(
    const zeroGradientAbsorbingWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::zeroGradientAbsorbingWallFvPatchScalarField::
zeroGradientAbsorbingWallFvPatchScalarField
(
    const zeroGradientAbsorbingWallFvPatchScalarField& wbppsf
)
:
    zeroGradientFvPatchScalarField(wbppsf)
{}


Foam::zeroGradientAbsorbingWallFvPatchScalarField::
zeroGradientAbsorbingWallFvPatchScalarField
(
    const zeroGradientAbsorbingWallFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroGradientAbsorbingWallFvPatchScalarField::updateCoeffs()
{
    DebugInFunction << endl;

    if (updated())
    {
        return;
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::snGrad() const
{
    DebugInFunction << endl;

    return this->patch().deltaCoeffs()*(scalar(0.0) - patchInternalField());
}

Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    DebugInFunction << endl;

    return tmp<scalarField>
    (
        new scalarField(this->size(), 1.0)
    );
}


Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    DebugInFunction << endl;

    return tmp<scalarField>
    (
        new scalarField(this->size(), Zero)
    );
}


Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::gradientInternalCoeffs()
const
{
    DebugInFunction << endl;

    return -this->patch().deltaCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::gradientBoundaryCoeffs()
const
{
    DebugInFunction << endl;

    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


void Foam::zeroGradientAbsorbingWallFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroGradientAbsorbingWallFvPatchScalarField
    );
}

// ************************************************************************* //
