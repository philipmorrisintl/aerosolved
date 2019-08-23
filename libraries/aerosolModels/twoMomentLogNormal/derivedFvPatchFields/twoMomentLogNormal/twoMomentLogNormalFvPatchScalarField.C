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

#include "twoMomentLogNormalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoMomentLogNormalFvPatchScalarField::
twoMomentLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    CMD_()
{}


Foam::twoMomentLogNormalFvPatchScalarField::
twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    CMD_(
        ptf.CMD_.valid()
        ?
        ptf.CMD_().clone().ptr()
        :
        NULL
    )
{}


Foam::twoMomentLogNormalFvPatchScalarField::
twoMomentLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false)
{
    CMD_ = Function1<scalar>::New("CMD", dict);

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
    else
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::twoMomentLogNormalFvPatchScalarField::
twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    CMD_(
        tppsf.CMD_.valid()
        ?
        tppsf.CMD_().clone().ptr()
        :
        NULL
    )
{}


Foam::twoMomentLogNormalFvPatchScalarField::
twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    CMD_(
        tppsf.CMD_.valid()
        ?
        tppsf.CMD_().clone().ptr()
        :
        NULL
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoMomentLogNormalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        const aerosolModel& aerosol =
            db().lookupObject<aerosolModel>("aerosolProperties");

        const aerosolThermo& thermo = aerosol.thermo();
        const speciesTable& dispSpecies = thermo.thermoDisp().species();

        const scalarField rhol(thermo.thermoDisp().rho(patch().index()));

        const scalar pi = constant::mathematical::pi;

        const scalar t = db().time().timeOutputValue();

        scalarField alpha(patch().size(), 0.0);

        forAll(dispSpecies, j)
        {
            alpha += thermo.Z()[j].boundaryField()[patch().index()];
        }

        const scalar sigma(aerosol.coeffs().lookupType<scalar>("sigma"));

        operator==
        (
            6.0*alpha*Foam::exp(-4.5*Foam::sqr(Foam::log(sigma)))
          / (pi*rhol*Foam::pow(CMD_->value(t),3.0))
        );
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::twoMomentLogNormalFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    CMD_->writeData(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        twoMomentLogNormalFvPatchScalarField
    );
}

// ************************************************************************* //
