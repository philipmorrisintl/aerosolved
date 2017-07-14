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

#include "twoMomentLogNormalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoMomentLogNormalFvPatchScalarField::twoMomentLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    sigma_(0),
    CMD_(0),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


Foam::twoMomentLogNormalFvPatchScalarField::twoMomentLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    sigma_(readScalar(dict.lookup("sigma"))),
    CMD_(readScalar(dict.lookup("CMD"))),
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


Foam::twoMomentLogNormalFvPatchScalarField::twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    sigma_(ptf.sigma_),
    CMD_(ptf.CMD_),
    phiName_(ptf.phiName_)
{}


Foam::twoMomentLogNormalFvPatchScalarField::twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    sigma_(tppsf.sigma_),
    CMD_(tppsf.CMD_),
    phiName_(tppsf.phiName_)
{}


Foam::twoMomentLogNormalFvPatchScalarField::twoMomentLogNormalFvPatchScalarField
(
    const twoMomentLogNormalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    sigma_(tppsf.sigma_),
    CMD_(tppsf.CMD_),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoMomentLogNormalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    checkName(name());

    const scalar pi = constant::mathematical::pi;

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

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

    if (aerosol.modType() != MOMENTAEROSOLMODEL)
    {
        FatalErrorIn
        (
            "Foam::twoMomentLogNormalFvPatchScalarField::updateCoeffs()"
        )   << "This boundary condition only works for a moment model." << nl
            << exit(FatalError);
    }

    const PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo.getProperty("rho", fluidThermo::LIQUID);

    // Compute liquid density and concentration

    scalarField ZtotRho(patch().size(), 0.0);
    scalarField Ztot(patch().size(), 0.0);
    scalarField rhoZ(patch().size(), 0.0);

    forAll(thermo.species(), j)
    {
        const word name(thermo.species().keys()[j] + "Z");

        const scalarField Z
        (
            max(patch().lookupPatchField<volScalarField, scalar>(name), 0.0)
        );

        forAll(ZtotRho, i)
        {
            ZtotRho[i] += Z[i] / stabilise(dataEntriesListRho_l[j].value(T[i]), SMALL);
        }

        Ztot += Z;

        rhoZ += rho*Z;
    }

    scalarField rhol(Ztot/stabilise(ZtotRho, SMALL));

    thermo.limitLiquidDensity(rhol);

    scalar d0(aerosol.params().lookupOrDefault<scalar>("Dcrit", 0.0));

    const scalar s(log(sigma_));

    const scalarField M
    (
        Ztot * 6.0/pi / rhol / pow(CMD_, 3.0) / exp(4.5*sqr(s))
      * (
            1.0 -
            (d0 > 0 ? erf((log(d0) - log(CMD_))/sqrt(2.0)/s) : -1)
        )
      / (
            1.0 -
            (d0 > 0 ? erf((log(d0) - log(CMD_) - 3.0*sqr(s))/sqrt(2.0)/s) : -1)
        )
    );

    this->refValue() = M;

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

inline Foam::word Foam::twoMomentLogNormalFvPatchScalarField::name()
{
    return this->dimensionedInternalField().name();
}

inline bool Foam::twoMomentLogNormalFvPatchScalarField::checkName
(
    const Foam::word name
)
{
    regExp r("M(_[0]+)?");

    if (!r.match(name))
    {
        FatalErrorIn("Foam::twoMomentLogNormalFvPatchScalarField::checkName()")
            << "This boundary conditions doens't work on a field named "
            << name << ". Only on M."
            << exit(FatalError);
    }

    return true;
}

void Foam::twoMomentLogNormalFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("sigma") << sigma_ << token::END_STATEMENT << nl;
    os.writeKeyword("CMD") << CMD_ << token::END_STATEMENT << nl;
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
