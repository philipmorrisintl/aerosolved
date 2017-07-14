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

#include "sectionalLogNormalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "aerosolModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalLogNormalFvPatchScalarField::sectionalLogNormalFvPatchScalarField
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


Foam::sectionalLogNormalFvPatchScalarField::sectionalLogNormalFvPatchScalarField
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


Foam::sectionalLogNormalFvPatchScalarField::sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& ptf,
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


Foam::sectionalLogNormalFvPatchScalarField::sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    sigma_(tppsf.sigma_),
    CMD_(tppsf.CMD_),
    phiName_(tppsf.phiName_)
{}


Foam::sectionalLogNormalFvPatchScalarField::sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    sigma_(tppsf.sigma_),
    CMD_(tppsf.CMD_),
    phiName_(tppsf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sectionalLogNormalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

    if (aerosol.modType() != SECTIONALAEROSOLMODEL)
    {
        FatalErrorIn
        (
            "Foam::sectionalLogNormalFvPatchScalarField::updateCoeffs()"
        )   << "This boundary condition only works for a sectional model." << nl
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

        scalarField Z
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

    const label ii(sectionNum());

    const label P(aerosol.P());

    const scalar s(log(sigma_));

    scalarField M(patch().size(), 0.0);

    forAll(M, k)
    {
        const scalarField d(pow(aerosol.x()*6.0/pi/rhol[k], 1.0/3.0));

        scalarField dz(P, 0.0);

        forAll(dz, i)
        {
            dz[i] = aerosol.y()[i+1]-aerosol.y()[i];
        }

        // Unscaled distribution transformed to mass space

        scalarField n
        (
            1.0/(sqrt(2.0*pi) * pow(d, 3.0) * s)
          * exp(-sqr(log(d)-log(CMD_))/(2.0*sqr(s)))
          * dz
        );

        // Scale distribution

        scalar Q = rhoZ[k] / sum(aerosol.x() * n);

        M[k] = n[ii] * Q / rho[k];
    }

    this->refValue() = M;

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

inline Foam::word Foam::sectionalLogNormalFvPatchScalarField::name()
{
    return this->dimensionedInternalField().name();
}

inline Foam::label Foam::sectionalLogNormalFvPatchScalarField::sectionNum()
{
    const word name(this->name());

    checkName(name);

    return readLabel(IStringStream(name.substr(2, name.find("_")))());
}

inline bool Foam::sectionalLogNormalFvPatchScalarField::checkName
(
    const Foam::word name
)
{
    regExp r("M\\.([0-9]+)(_[0]+)?");

    if (!r.match(name))
    {
        FatalErrorIn("Foam::sectionalLogNormalFvPatchScalarField::checkName()")
            << "This boundary conditions doens't work on a field named "
            << name
            << exit(FatalError);
    }

    return true;
}

void Foam::sectionalLogNormalFvPatchScalarField::write(Ostream& os) const
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
        sectionalLogNormalFvPatchScalarField
    );
}

// ************************************************************************* //
