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

#include "saturatedMixtureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::saturatedMixtureFvPatchScalarField::checkFields()
{
    if (!checkFields_)
    {
        const aerosolThermo& thermo =
            db().lookupObject<aerosolThermo>("thermophysicalProperties");

        forAll(thermo.contSpecies(), j)
        {
            if
            (
                thermo.Y()[j].boundaryField()[patch().index()].type()
             != type()
            )
            {
                FatalErrorInFunction
                    << "All continuous species should have the "
                    << "saturatedMixture BC at this patch" << nl
                    << exit(FatalError);
            }
        }

        checkFields_ = true;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    S_(ptf.S_),
    f_(ptf.f_)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false)
{
    if (dict.found("S"))
    {
        const List<keyType> SKeys(dict.subDict("S").keys());

        forAll(SKeys, i)
        {
            S_.insert
            (
                SKeys[i],
                readScalar(dict.subDict("S").lookup(SKeys[i]))
            );
        }
    }

    const List<keyType> fKeys(dict.subDict("inertMoleFrac").keys());

    forAll(fKeys, i)
    {
        f_.insert
        (
            fKeys[i],
            readScalar(dict.subDict("inertMoleFrac").lookup(fKeys[i]))
        );
    }

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


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    S_(tppsf.S_),
    f_(tppsf.f_)
{}


Foam::saturatedMixtureFvPatchScalarField::
saturatedMixtureFvPatchScalarField
(
    const saturatedMixtureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    S_(tppsf.S_),
    f_(tppsf.f_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::saturatedMixtureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    checkFields();

    if (db().foundObject<aerosolThermo>("thermophysicalProperties"))
    {
        aerosolThermo& thermo =
            db().lookupObjectRef<aerosolThermo>("thermophysicalProperties");

        rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

        const speciesTable& contSpecies = thermo.contSpecies();
        const speciesTable& activeSpecies = thermo.activeSpecies();
        const speciesTable& inactiveSpecies = thermo.inactiveSpecies();

        const scalarField& p = thermo.p().boundaryField()[patch().index()];
        const scalarField& T = thermo.T().boundaryField()[patch().index()];

        scalarField Y(patch().size(), 0.0);

        const word specieName(internalField().member());

        const label jj(contSpecies[specieName]);

        scalarList W(contSpecies.size(), 0.0);

        forAll(contSpecies, j)
        {
            W[j] = thermoCont.composition().W(j);
        }

        PtrList<scalarField> pSat(activeSpecies.size());

        scalarList S(activeSpecies.size(), 1.0);

        forAll(activeSpecies, j)
        {
            pSat.set
            (
                j,
                new scalarField
                (
                    thermoCont.property(activeSpecies[j], "pSat").value(T)
                )
            );

            if (S_.found(activeSpecies[j]))
            {
                S[j] = S_[activeSpecies[j]];
            }
        }

        scalarList inertMoleFrac(inactiveSpecies.size(), 0.0);

        forAll(inactiveSpecies, j)
        {
            inertMoleFrac[j] = f_[inactiveSpecies[j]];
        }

        inertMoleFrac = inertMoleFrac/sum(inertMoleFrac);

        forAll(*this, facei)
        {
            scalarList X
            (
                S
              * entryList(pSat,facei)
              / p[facei]
            );

            const scalar sumX(min(sum(X),1.0));

            forAll(inactiveSpecies, j)
            {
                X.append((1.0-sumX)*inertMoleFrac[j]);
            }

            X = X/sum(X);

            Y[facei] = X[jj]*W[jj]/sum(X*W);
        }

        operator==(Y);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::saturatedMixtureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);

    if (S_.size() > 0)
    {
        os.beginBlock("S");

        forAll(S_.toc(), i)
        {
            const word key(S_.toc()[i]);

            os.writeKeyword(key)
                << token::SPACE << S_[key] << token::END_STATEMENT << nl;
        }

        os.endBlock();
    }

    os.beginBlock("inertMoleFrac");

    forAll(f_.toc(), i)
    {
        const word key(f_.toc()[i]);

        os.writeKeyword(key)
            << token::SPACE << f_[key] << token::END_STATEMENT << nl;
    }

    os.endBlock();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        saturatedMixtureFvPatchScalarField
    );
}

// ************************************************************************* //
