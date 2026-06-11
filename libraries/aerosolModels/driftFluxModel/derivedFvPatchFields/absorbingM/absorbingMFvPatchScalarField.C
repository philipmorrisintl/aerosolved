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

#include <iostream>
#include "absorbingMFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "driftFluxModel.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absorbingMFvPatchScalarField::
absorbingMFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    alpha_(1.0),
    beta_(1.0),
    clipInwardFlux_(true)
{}


Foam::absorbingMFvPatchScalarField::
absorbingMFvPatchScalarField
(
    const absorbingMFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    alpha_(ptf.alpha_),
    beta_(ptf.beta_),
    clipInwardFlux_(ptf.clipInwardFlux_)
{}


Foam::absorbingMFvPatchScalarField::
absorbingMFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    alpha_(dict.lookupOrDefault<scalar>("reboundValue", 1.0)),
    beta_(dict.lookupOrDefault<scalar>("stickinessValue", 1.0)),
    clipInwardFlux_(dict.lookupOrDefault<bool>("clipInwardFlux", true))
{}


Foam::absorbingMFvPatchScalarField::
absorbingMFvPatchScalarField
(
    const absorbingMFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    alpha_(tppsf.alpha_),
    beta_(tppsf.beta_),
    clipInwardFlux_(tppsf.clipInwardFlux_)
{}


Foam::absorbingMFvPatchScalarField::
absorbingMFvPatchScalarField
(
    const absorbingMFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    alpha_(tppsf.alpha_),
    beta_(tppsf.beta_),
    clipInwardFlux_(tppsf.clipInwardFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::label Foam::absorbingMFvPatchScalarField::sectionIndex()
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
        const word sectioNumi(system.distribution()[i].M().group());

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
void Foam::absorbingMFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().foundObject<aerosolModel>("aerosolProperties"))
    {
        const aerosolModel& aerosol =
            db().lookupObject<aerosolModel>("aerosolProperties");

        const scalarField phiUd
        (
            aerosol.phi().boundaryField()[patch().index()]
          - aerosol.drift().phiCorr().boundaryField()[patch().index()]
        );

        const scalarField rho(aerosol.rho().boundaryField()[patch().index()]);

        // Get inertial drift flux and Brownian diffusivity from the sectional
        // system, or if that does not exist from the aerosol system

        scalarField phiInertial(patch().size());
        scalarField D(patch().size());

        if (db().foundObject<fixedSectionalSystem>("fixedSectionalSystem"))
        {
            const fixedSectionalSystem& system =
                db().lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

            const label sec(sectionIndex());

            phiInertial =
                system.distribution()[sec].phiInertial()
               .boundaryField()[patch().index()];

            D = system.distribution()[sec].D()
               .boundaryField()[patch().index()];
        }
        else
        {
            phiInertial =
                aerosol.drift().phiInertial()
               .boundaryField()[patch().index()];

            D = aerosol.drift().DDisp()
               .boundaryField()[patch().index()];
        }

        // Dispersed velocity

        const scalarField Ud(phiUd/(rho*patch().magSf()));

        // Brownian velocity

        const scalarField VBrownian(-D/(patch().delta() & patch().nf()));

        // Inertial velocity

        const scalarField VInertial(phiInertial/(rho*patch().magSf()));

        // Get new patch coefficients

        const scalarField Mf
        (
            patchInternalField()
          * aerosol.rho().boundaryField()[patch().index()].patchInternalField()
          / rho
          * (
                alpha_
              - (alpha_ + beta_ - scalar(1.0))
              / (scalar(1.0) - mag(Ud+VInertial)/(VBrownian+SMALL))
            )
        );

        // Check for Inward Flux

        if (clipInwardFlux_)
        {
            operator==(Mf*pos(Ud+VInertial));
        }
        else
        {
            operator==(Mf);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::absorbingMFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntryIfDifferent<scalar>("reboundValue", 1.0, alpha_);
    os.writeEntryIfDifferent<scalar>("stickinessValue", 1.0, beta_);
    os.writeEntryIfDifferent<bool>("clipInwardFlux", true, clipInwardFlux_);

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        absorbingMFvPatchScalarField
    );
}

// ************************************************************************* //
