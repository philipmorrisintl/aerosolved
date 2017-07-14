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

#include "saturatedVaporMixtureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    checkFields_(false),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
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

    const Switch phaseChange = thermo.phaseChange()[j];

    if (!phaseChange)
    {
        const word speciesName = thermo.species().keys()[j];

        if (!thermo.species().subDict(speciesName).found("inertFraction"))
        {
            FatalErrorIn
            (
                "Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField(...)"
            )   << "inertFraction keyword not found for " << speciesName << ". "
                << "Required for this boundary condition." << nl
                << exit(FatalError);
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField
(
    const saturatedVaporMixtureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    checkFields_(false),
    phiName_(ptf.phiName_)
{}


Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField
(
    const saturatedVaporMixtureFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    checkFields_(false),
    phiName_(tppsf.phiName_)
{}


Foam::saturatedVaporMixtureFvPatchScalarField::saturatedVaporMixtureFvPatchScalarField
(
    const saturatedVaporMixtureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    checkFields_(false),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::saturatedVaporMixtureFvPatchScalarField::updateCoeffs()
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

    // Find the necessary references

    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    const fvPatchField<scalar>& p1 =
        patch().lookupPatchField<volScalarField, scalar>("p1");

    PtrList<DataEntry<scalar> >& dataEntriesListP_s =
        thermo.getProperty("P_s", fluidThermo::VAPOR);

    const scalarList& M = thermo.M();

    const scalar p0 = thermo.p0().value();

    // Compute the mole fractions for all species

    List< List<scalar> > X
    (
        thermo.nSpecies(),
        List<scalar>(patch().size(), 0.0)
    );

    scalarField XphaseChangeSum(patch().size(), 0.0);

    // Phase changing species

    forAll(thermo.speciesPhaseChange(), j)
    {
        // Set the equilibrium vapor pressure for each face

        scalarField pjeq(patch().size(), 0.0);

        forAll(*this, faceI)
        {
            pjeq[faceI] = dataEntriesListP_s[j].value(T[faceI]);
        }

        // Saturation mole fraction

        X[j] = pjeq / (p0 + p1);

        XphaseChangeSum += X[j];
    }

    // Inert species

    forAll(thermo.species(), j)
    {
        if (!thermo.phaseChange()[j])
        {
            const word speciesName = thermo.species().keys()[j];

            scalar inertFraction(0.0);

            thermo.species().subDict(speciesName).lookup("inertFraction") >>
                inertFraction;

            X[j] = inertFraction * (1.0 - XphaseChangeSum);
        }
    }

    scalarField Mmix(patch().size(), 0.0);

    forAll(thermo.species(), j)
    {
        Mmix += X[j] * M[j];
    }

    const label j = speciesIndex();

    Yj = X[j] * M[j] / Mmix;

    this->refValue() = Yj;

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

void Foam::saturatedVaporMixtureFvPatchScalarField::checkFields()
{
    if (!checkFields_)
    {
        fluidThermo& thermo =
            const_cast<fluidThermo&>
            (
                db().lookupObject<fluidThermo>("fluidThermoProperties")
            );

        if (thermo.nSpecies() - thermo.nSpeciesPhaseChange() == 0)
        {
            FatalErrorIn
            (
                "Foam::saturatedVaporMixtureFvPatchScalarField::checkFields()"
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
                    "Foam::saturatedVaporMixtureFvPatchScalarField::checkFields()"
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
                    "Foam::saturatedVaporMixtureFvPatchScalarField::checkFields()"
                )   << "The sum of the inertFraction values for all non-phase "
                    << "changing species should be equal to one but is instead "
                    << "equal to " << f << "." << exit(FatalError);
        }

        checkFields_ = true;
    }
}

Foam::label Foam::saturatedVaporMixtureFvPatchScalarField::speciesIndex()
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
            "Foam::saturatedVaporMixtureFvPatchScalarField::updateCoeffs()"
        )   << "Could not find species index number for field " << name << nl
            << exit(FatalError);
    }

    return j;
}

void Foam::saturatedVaporMixtureFvPatchScalarField::write(Ostream& os) const
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
        saturatedVaporMixtureFvPatchScalarField
    );
}

// ************************************************************************* //
