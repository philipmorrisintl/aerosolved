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

#include "zeroGradientAbsorbingWallFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

Foam::tmp<Foam::scalarField>
Foam::zeroGradientAbsorbingWallFvPatchScalarField::snGrad() const
{
    return this->patch().deltaCoeffs()*(scalar(0.0) - patchInternalField());
}

void Foam::zeroGradientAbsorbingWallFvPatchScalarField::write(Ostream& os) const
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
