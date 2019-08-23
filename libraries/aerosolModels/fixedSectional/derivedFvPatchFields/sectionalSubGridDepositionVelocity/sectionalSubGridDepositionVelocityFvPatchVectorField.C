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

#include "sectionalSubGridDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "aerosolModel.H"
#include "subGridDepositionModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::sectionIndex()
const
{
    const fixedSectionalSystem& system =
        db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    const word sectionNum
    (
        this->internalField().group()
    );

    forAll(system.distribution(), i)
    {
        const word sectioNumi
        (
            system.distribution()[i].M().group()
        );

        if (sectionNum == sectioNumi)
        {
            return i;
        }
    }

    FatalErrorInFunction
        << "The section to which this BC's field ("
        << this->internalField().name()
        << ") belongs could not be found"
        << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    g_(vector::zero),
    maxIter_(99),
    tolerance_(1e-8)
{}


Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const sectionalSubGridDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    g_(ptf.g_),
    maxIter_(ptf.maxIter_),
    tolerance_(ptf.tolerance_)
{}


Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF),
    g_
    (
        uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                patch().boundaryMesh().mesh().time().constant(),
                patch().boundaryMesh().mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedVector("g", dimVelocity/dimTime, vector::zero)
        ).value()
    ),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 99)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-8))
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::
sectionalSubGridDepositionVelocityFvPatchVectorField
(
    const sectionalSubGridDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(fcvpvf, iF),
    g_(fcvpvf.g_),
    maxIter_(fcvpvf.maxIter_),
    tolerance_(fcvpvf.tolerance_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    const fixedSectionalSystem& system =
        db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    const aerosolModel& aerosol = system.aerosol();

    const section& sec = system.distribution()[sectionIndex()];

    const vectorField Up
    (
        aerosol.U().boundaryField()[patch().index()].patchInternalField()
    );

    const aerosolThermo& thermo = aerosol.thermo();

    const scalarField rhoc(thermo.thermoCont().rho(patch().index()));
    const scalarField rhod(thermo.thermoDisp().rho(patch().index()));

    const scalarField gamma(rhoc/rhod);

    const scalarField mu(thermo.thermoCont().mu(patch().index()));

    const scalarField d(sec.d(rhod));

    const scalarField tau(rhod*sqr(d)/(18.0*mu));

    const scalarField delta(mag(patch().delta()));
    const vectorField n(patch().nf());

    const scalarField v0(-((patchInternalField() + Up) & n) * tau/delta);

    const scalarField g(-(g_ & n) * (1.0-gamma)*sqr(tau)/delta);

    const scalarField u(-(Up & n) * tau/delta);

    scalarField v(patch().size(), 0.0);

    forAll(*this, facei)
    {
        if (v0[facei] < tolerance_)
        {
            if (delta[facei] > d[facei]/2.0)
            {
                const subGridDepositionModel model
                (
                    u[facei],
                    g[facei],
                    v0[facei],
                    d[facei]/2.0,
                    maxIter_,
                    tolerance_
                );

                if (model.collision())
                {
                    const scalar t(model.t());
                    v[facei] = model.v(t);
                }
            }
            else
            {
                v[facei] = v0[facei];
            }
        }
    }

    operator==(-v*delta/tau*n);
}


void Foam::sectionalSubGridDepositionVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);

    os.writeEntryIfDifferent<label>("maxIter", 99, maxIter_);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-8, tolerance_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sectionalSubGridDepositionVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
