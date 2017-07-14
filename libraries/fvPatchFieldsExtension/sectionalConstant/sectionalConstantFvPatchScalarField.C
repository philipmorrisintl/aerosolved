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

#include "sectionalConstantFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalConstantFvPatchScalarField::sectionalConstantFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::sectionalConstantFvPatchScalarField::sectionalConstantFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    this->refValue() = pTraits<scalar>::zero;

    // Read rhol, for later use

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    thermo.readProperty("rho", fluidThermo::LIQUID, thermo.speciesPhaseChange());

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::sectionalConstantFvPatchScalarField::sectionalConstantFvPatchScalarField
(
    const sectionalConstantFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


Foam::sectionalConstantFvPatchScalarField::sectionalConstantFvPatchScalarField
(
    const sectionalConstantFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    phiName_(tppsf.phiName_)
{}


Foam::sectionalConstantFvPatchScalarField::sectionalConstantFvPatchScalarField
(
    const sectionalConstantFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sectionalConstantFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    // Get a reference to the thermo and aerosol objects

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    aerosolModel& aerosol =
        const_cast<aerosolModel&>
        (
            db().lookupObject<aerosolModel>("aerosolProperties")
        );

    if (aerosol.modType() != SECTIONALAEROSOLMODEL)
    {
        FatalErrorIn
        (
            "Foam::sectionalLogNormalFvPatchScalarField::updateCoeffs()"
        )   << "This boundary condition only works for a sectional model." << nl
            << exit(FatalError);
    }

    // Compute total Z

    scalarField Ztot(patch().size(), 0.0);

    forAll(thermo.species(), j)
    {
        const word name(thermo.species().keys()[j] + "Z");

        const fvPatchField<scalar>& Z =
            patch().lookupPatchField<volScalarField, scalar>(name);

        Ztot += Z;
    }

    this->refValue() = Ztot/sum(aerosol.x());

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

inline Foam::word Foam::sectionalConstantFvPatchScalarField::name()
{
    return this->dimensionedInternalField().name();
}

inline Foam::label Foam::sectionalConstantFvPatchScalarField::sectionNum()
{
    const word name(this->name());

    checkName(name);

    return readLabel(IStringStream(name.substr(2, name.find("_")))());
}

inline bool Foam::sectionalConstantFvPatchScalarField::checkName
(
    const Foam::word name
)
{
    regExp r("M\\.([0-9]+)(_[0]+)?");

    if (!r.match(name))
    {
        FatalErrorIn("Foam::sectionalConstantFvPatchScalarField::checkName()")
            << "This boundary conditions doens't work on a field named "
            << name
            << exit(FatalError);
    }

    return true;
}

void Foam::sectionalConstantFvPatchScalarField::write(Ostream& os) const
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
        sectionalConstantFvPatchScalarField
    );
}

// ************************************************************************* //
