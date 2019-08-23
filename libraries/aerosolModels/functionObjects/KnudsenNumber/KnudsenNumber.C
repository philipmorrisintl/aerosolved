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

#include "KnudsenNumber.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(KnudsenNumber, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        KnudsenNumber,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::KnudsenNumber::calc()
{
    const tmp<volScalarField> tKnudsenNumber
    (
        new volScalarField
        (
            IOobject
            (
                "t" + resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            2.0*aerosol_.thermo().lambda()/aerosol_.meanDiameter(p_,q_)
        )
    );

    return store(resultName_, tKnudsenNumber, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::KnudsenNumber::KnudsenNumber
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    p_(dict.lookupOrDefault<scalar>("p", 1)),
    q_(dict.lookupOrDefault<scalar>("q", 0)),
    aerosol_(lookupObject<aerosolModel>("aerosolProperties")),
    resultName_
    (
        dict.lookupOrDefault<word>
        (
            "result",
            word("KnudsenNumber("+Foam::name(p_)+","+Foam::name(q_)+")")
        )
    )
{
    calc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::KnudsenNumber::~KnudsenNumber()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::KnudsenNumber::execute()
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


bool Foam::functionObjects::KnudsenNumber::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::KnudsenNumber::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
