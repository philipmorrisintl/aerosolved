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

#include "fixedMeanVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedMeanVelocityFvPatchVectorField::
fixedMeanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    meanU_(),
    curTimeIndex_(-1)
{}


Foam::fixedMeanVelocityFvPatchVectorField::
fixedMeanVelocityFvPatchVectorField
(
    const fixedMeanVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    meanU_(ptf.meanU_().clone().ptr()),
    curTimeIndex_(-1)
{}


Foam::fixedMeanVelocityFvPatchVectorField::
fixedMeanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    curTimeIndex_(-1)
{
    meanU_ = DataEntry<scalar>::New("meanU", dict);

    evaluate(Pstream::blocking);
}


Foam::fixedMeanVelocityFvPatchVectorField::
fixedMeanVelocityFvPatchVectorField
(
    const fixedMeanVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    meanU_(ptf.meanU_().clone().ptr()),
    curTimeIndex_(-1)
{}


Foam::fixedMeanVelocityFvPatchVectorField::
fixedMeanVelocityFvPatchVectorField
(
    const fixedMeanVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    meanU_(ptf.meanU_().clone().ptr()),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedMeanVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().value();

    const scalar meanU(meanU_->value(t));

    const vectorField Ui(patchInternalField());

    const scalar meanUi = gSum(Ui & patch().Sf()) / gSum(patch().magSf());

    vectorField Up(Ui);

    if (meanU <= 0.0)
    {
        // Push plug flow

        Up = meanU * patch().nf();
    }
    else
    {
        if (mag(meanU) > SMALL && mag(meanUi)/mag(meanU) > 0.5)
        {
            // Pull the internal field profile

            Up = Up/meanUi * meanU;
        }
        else
        {
            // Plug flow

            Up = meanU * patch().nf();
        }

    }

    curTimeIndex_ = this->db().time().timeIndex();

    operator==(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedMeanVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
    meanU_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       fixedMeanVelocityFvPatchVectorField
   );
}

// ************************************************************************* //
