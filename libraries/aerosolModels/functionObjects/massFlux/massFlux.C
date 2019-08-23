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

#include "massFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        massFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::massFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Mass fluxes [kg/s]");
    writeCommented(os, "Time");

    forAll(fluxFieldNames_, fluxi)
    {
        const word fluxFieldName(fluxFieldNames_[fluxi]);

        writeTabbed(os, fluxFieldName);
    }

    os  << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massFlux::massFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sampleFlux(name, runTime, dict)
{
    const aerosolThermo& thermo =
        mesh_.lookupObject<aerosolThermo>("thermophysicalProperties");

    fluxFieldNames_.setSize(thermo.Y().size()+thermo.Z().size());

    forAll(thermo.Y(), j)
    {
        fluxFieldNames_[j] =
            IOobject::groupName("phiEff", thermo.Y()[j].name());
    }

    forAll(thermo.Z(), j)
    {
        const label k(j+thermo.Y().size());

        fluxFieldNames_[k] =
            IOobject::groupName("phiEff", thermo.Z()[j].name());
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massFlux::~massFlux()
{}

// ************************************************************************* //
