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

#include "sectionalLogNormalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sectionalLogNormalFvPatchScalarField::sectionIndex() const
{
    return readLabel(IStringStream(internalField().group())());
}

Foam::scalarField Foam::sectionalLogNormalFvPatchScalarField::logNormalIntegral
(
    const scalar& xl,
    const scalar& xu,
    const scalarField& CMM
) const
{
    return 0.5*
    (
        erf(log(xu/CMM)/(3.0*sqrt(2.0)*log(sigma_)))
      - erf(log(max(xl/CMM,VSMALL))/(3.0*sqrt(2.0)*log(sigma_)))
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalLogNormalFvPatchScalarField::
sectionalLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    sigma_(0.0),
    gamma_(0.0),
    CMD_()
{}


Foam::sectionalLogNormalFvPatchScalarField::
sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    sigma_(ptf.sigma_),
    gamma_(ptf.gamma_),
    CMD_(
        ptf.CMD_.valid()
        ?
        ptf.CMD_().clone().ptr()
        :
        NULL
    )
{}


Foam::sectionalLogNormalFvPatchScalarField::
sectionalLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    sigma_(readScalar(dict.lookup("sigma"))),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 3.0)),
    CMD_()
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


Foam::sectionalLogNormalFvPatchScalarField::
sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    sigma_(tppsf.sigma_),
    gamma_(tppsf.gamma_),
    CMD_(
        tppsf.CMD_.valid()
        ?
        tppsf.CMD_().clone().ptr()
        :
        NULL
    )
{}


Foam::sectionalLogNormalFvPatchScalarField::
sectionalLogNormalFvPatchScalarField
(
    const sectionalLogNormalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    sigma_(tppsf.sigma_),
    gamma_(tppsf.gamma_),
    CMD_(
        tppsf.CMD_.valid()
        ?
        tppsf.CMD_().clone().ptr()
        :
        NULL
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sectionalLogNormalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().foundObject<fixedSectionalSystem>("fixedSectionalSystem"))
    {
        const scalar pi = constant::mathematical::pi;

        const fixedSectionalSystem& system =
            db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

        const aerosolModel& aerosol = system.aerosol();
        const aerosolThermo& thermo = aerosol.thermo();
        const sectionalDistribution& dist = system.distribution();
        const speciesTable& dispSpecies = thermo.thermoDisp().species();

        const scalar t = db().time().timeOutputValue();

        scalarField alpha(patch().size(), 0.0);

        forAll(dispSpecies, j)
        {
            alpha += thermo.Z()[j].boundaryField()[patch().index()];
        }

        const scalarField rhol(thermo.thermoDisp().rho(patch().index()));

        const scalarField CMM(pi/6.0*rhol*Foam::pow(CMD_->value(t), 3.0));

        scalarField F(patch().size(), 0.0);

        forAll(dist, i)
        {
            const scalarField di(pow(6.0*dist.x()[i]/(pi*rhol), 1.0/3.0));

            const scalarField Fi
            (
                logNormalIntegral(dist.y()[i], dist.y()[i+1], CMM)
            );

            F += pow(di,gamma_)*Fi;
        }

        const scalarField a
        (
            pow(CMD_->value(t),gamma_-3.0)/F
          * 6.0/pi*alpha/rhol
          * exp((0.5*sqr(gamma_)-4.5)*sqr(log(sigma_)))
        );

        const label i(sectionIndex()-1);

        const scalarField Fi
        (
            logNormalIntegral(dist.y()[i], dist.y()[i+1], CMM)
        );

        operator==(a*Fi);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::sectionalLogNormalFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;
    CMD_->writeData(os);
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
