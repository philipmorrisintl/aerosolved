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

#include "thermoField.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermoField, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        thermoField,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::thermoField::thermoFieldType
>
Foam::functionObjects::thermoField::thermoFieldTypeNames_
{
    { thermoFieldType::Cp, "Cp" },
    { thermoFieldType::Cv, "Cv" },
    { thermoFieldType::Cpv, "Cpv" },
    { thermoFieldType::CpByCpv, "CpByCpv" },
    { thermoFieldType::gamma, "gamma" },
    { thermoFieldType::hc, "hc" },
    { thermoFieldType::nu, "nu" },
    { thermoFieldType::kappa, "kappa" },
    { thermoFieldType::sumY, "sumY" },
    { thermoFieldType::sumZ, "sumZ" },
    { thermoFieldType::lambda, "lambda" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::thermoField::calc()
{
    IOobject fHeader
    (
        "t" + resultName_,
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    switch (thermoField_)
    {
        case Cp:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.Cp())
            );

            return store(resultName_, tf, true);
        }

        case Cv:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.Cv())
            );

            return store(resultName_, tf, true);
        }

        case Cpv:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.Cpv())
            );

            return store(resultName_, tf, true);
        }

        case CpByCpv:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.CpByCpv())
            );

            return store(resultName_, tf, true);
        }

        case gamma:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.gamma())
            );

            return store(resultName_, tf, true);
        }

        case hc:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.hc())
            );

            return store(resultName_, tf, true);
        }

        case nu:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.nu())
            );

            return store(resultName_, tf, true);
        }

        case kappa:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.kappa())
            );

            return store(resultName_, tf, true);
        }

        case sumY:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.sumY())
            );

            return store(resultName_, tf, true);
        }

        case sumZ:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.sumZ())
            );

            return store(resultName_, tf, true);
        }

        case lambda:
        {
            const tmp<volScalarField> tf
            (
                new volScalarField(fHeader, thermo_.lambda())
            );

            return store(resultName_, tf, true);
        }

        default:

            return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermoField::thermoField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    thermoField_(thermoFieldTypeNames_.lookup("thermoField", dict)),
    thermo_(lookupObjectRef<aerosolThermo>("thermophysicalProperties")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            thermoFieldTypeNames_[thermoField_]
        )
    )
{
    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::thermoField::~thermoField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::thermoField::execute()
{
    if (!calc())
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        // Clear the result field from the objectRegistry if present
        clear();

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::thermoField::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::thermoField::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
