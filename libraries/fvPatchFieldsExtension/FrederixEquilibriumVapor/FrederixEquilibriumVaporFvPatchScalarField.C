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

#include "FrederixEquilibriumVaporFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    phaseChange_(false),
    inertFraction_(0.0),
    checkFields_(false),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    phaseChange_(false),
    inertFraction_(0.0),
    checkFields_(false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    this->refValue() = pTraits<scalar>::zero;

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    thermo.readProperty("P_s", fluidThermo::VAPOR, thermo.speciesPhaseChange());

    const label j = speciesIndex();

    phaseChange_ = thermo.phaseChange()[j];

    if (!phaseChange_)
    {
        const word speciesName = thermo.species().keys()[j];

        if (!thermo.species().subDict(speciesName).found("inertFraction"))
        {
            FatalErrorIn
            (
                "Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField(...)"
            )   << "inertFraction keyword not found for " << speciesName << ". "
                << "Required for this boundary condition." << nl
                << exit(FatalError);
        }

        thermo.species().subDict(speciesName).lookup("inertFraction") >>
            inertFraction_;
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField
(
    const FrederixEquilibriumVaporFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    phaseChange_(ptf.phaseChange_),
    inertFraction_(ptf.inertFraction_),
    checkFields_(ptf.checkFields_),
    phiName_(ptf.phiName_)
{}


Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField
(
    const FrederixEquilibriumVaporFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    phaseChange_(tppsf.phaseChange_),
    inertFraction_(tppsf.inertFraction_),
    checkFields_(tppsf.checkFields_),
    phiName_(tppsf.phiName_)
{}


Foam::FrederixEquilibriumVaporFvPatchScalarField::FrederixEquilibriumVaporFvPatchScalarField
(
    const FrederixEquilibriumVaporFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    phaseChange_(tppsf.phaseChange_),
    inertFraction_(tppsf.inertFraction_),
    checkFields_(tppsf.checkFields_),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FrederixEquilibriumVaporFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    checkFields();

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalarField Yj(patch().size(), 0.0);

    if (phaseChange_)
    {
        // Find the necessary references

        const fvPatchField<scalar>& T =
            patch().lookupPatchField<volScalarField, scalar>("T");

        const fvPatchField<scalar>& p1 =
            patch().lookupPatchField<volScalarField, scalar>("p1");

        PtrList<DataEntry<scalar> >& dataEntriesListP_s =
            thermo.getProperty("P_s", fluidThermo::VAPOR);

        const scalarList& M = thermo.M();

        const scalar p0 = thermo.p0().value();

        const label j = speciesIndex();

        const word speciesName = thermo.species().keys()[j];

        // References for this species

        const scalarField& Zj =
            patch().lookupPatchField<volScalarField, scalar>(word(speciesName + "Z"));

        const scalar Mj = M[j];

        // Compute the j mole fraction w.r.t. liquid phase

        scalarField sumZkOverMk(patch().size(), 0.0);

        forAll(thermo.species(), k)
        {
            const word name(thermo.species().keys()[k] + "Z");

            const scalarField& Zk =
                patch().lookupPatchField<volScalarField, scalar>(name);

            sumZkOverMk += max(Zk, 0.0) / M[k];
        }

        const scalarField omegaj = max(Zj, 0.0)/Mj/stabilise(sumZkOverMk, SMALL);

        // Set the equilibrium vapor pressure for each face

        scalarField pjeq(patch().size(), 0.0);

        forAll(*this, faceI)
        {
            pjeq[faceI] = dataEntriesListP_s[j].value(T[faceI]);
        }

        // Compute the hypothetical equilibrium mole fraction

        const scalarField chijeq = omegaj * pjeq / (p0 + p1);

        // Compute the vapor mean molecular wait

        scalarField sumYkOverMk(patch().size(), 0.0);
        scalarField sumYk(patch().size(), 0.0);

        forAll(thermo.species(), k)
        {
            const word name(thermo.species().keys()[k] + "Y");

            const scalarField& Yk =
                patch().lookupPatchField<volScalarField, scalar>(name);

            sumYkOverMk += max(Yk, 0.0) / M[k];

            sumYk += max(Yk, 0.0);
        }

        const scalarField Mv = sumYk / stabilise(sumYkOverMk, SMALL);

        // Set the equilibrium vapor mass fraction

        scalarField Mjeq = (1.0 - chijeq)*Mv + chijeq*Mj;

        Yj = chijeq * Mj / Mjeq;
    }
    else
    {
        Yj = 1.0;

        forAll(thermo.speciesPhaseChange(), k)
        {
            const word name(thermo.species().keys()[k] + "Y");

            const scalarField& Yk =
                patch().lookupPatchField<volScalarField, scalar>(name);

            Yj -= Yk;
        }

        forAll(thermo.species(), k)
        {
            const word name(thermo.species().keys()[k] + "Z");

            const scalarField& Zk =
                patch().lookupPatchField<volScalarField, scalar>(name);

            Yj -= Zk;
        }

        Yj *= inertFraction_;
    }

    this->refValue() = Yj;

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

void Foam::FrederixEquilibriumVaporFvPatchScalarField::checkFields()
{
    if (!checkFields_)
    {
        fluidThermo& thermo =
            const_cast<fluidThermo&>
            (
                db().lookupObject<fluidThermo>("fluidThermoProperties")
            );

        if (thermo.nSpeciesPhaseChange() == 0)
        {
            FatalErrorIn
            (
                "Foam::FrederixEquilibriumVaporFvPatchScalarField::checkFields()"
            )   << "No phase changing species found. This boundary condition "
                << "needs at least one." << nl << exit(FatalError);
        }

        if (thermo.nSpecies() - thermo.nSpeciesPhaseChange() == 0)
        {
            FatalErrorIn
            (
                "Foam::FrederixEquilibriumVaporFvPatchScalarField::checkFields()"
            )   << "No non-phase changing species found. This boundary condition "
                << "needs at least one." << nl << exit(FatalError);
        }

        scalar f(0.0);
        scalar fk(0.0);

        forAll(thermo.species(), k)
        {
            const word name(thermo.species().keys()[k] + "Y");
            const word speciesName(name.substr(0, name.length()-1));

            const fvPatchScalarField& Yk =
                patch().lookupPatchField<volScalarField, scalar>(name);

            if (Yk.type() != type())
            {
                FatalErrorIn
                (
                    "Foam::FrederixEquilibriumVaporFvPatchScalarField::checkFields()"
                )   << "Not all Y fields have the " << type()
                    << " boundary condition on patch "
                    << patch().name() << nl << exit(FatalError);
            }

            if (!thermo.phaseChange()[k])
            {
                thermo.species().subDict(speciesName).lookup("inertFraction") >> fk;
                f += fk;
            }
        }

        if (f < 1.0-SMALL || f > 1.0+SMALL)
        {
                FatalErrorIn
                (
                    "Foam::FrederixEquilibriumVaporFvPatchScalarField::checkFields()"
                )   << "The sum of the inertFraction values for all non-phase "
                    << "changing species should be equal to one but is instead "
                    << "equal to " << f << "." << exit(FatalError);
        }

        checkFields_ = true;
    }
}

Foam::label Foam::FrederixEquilibriumVaporFvPatchScalarField::speciesIndex()
{
    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    const word name(this->dimensionedInternalField().name());
    const word speciesName(name.substr(0, name.length()-1));

    label j(-1);

    forAll(thermo.species(), jj)
    {
        if (thermo.species().keys()[jj] == speciesName)
        {
            j = jj;
            break;
        }
    }

    if (j == -1)
    {
        FatalErrorIn
        (
            "Foam::FrederixEquilibriumVaporFvPatchScalarField::updateCoeffs()"
        )   << "Could not find species index number for field " << name << nl
            << exit(FatalError);
    }

    return j;
}

void Foam::FrederixEquilibriumVaporFvPatchScalarField::write(Ostream& os) const
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
        FrederixEquilibriumVaporFvPatchScalarField
    );
}

// ************************************************************************* //
