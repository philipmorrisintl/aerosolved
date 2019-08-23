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

#include "sectionalFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sectionalFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sectionalFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sectionalFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Sectional number flux [#/s]");
    writeCommented(os, "Time");

    forAll(fluxFieldNames_, fluxi)
    {
        const word fluxFieldName(fluxFieldNames_[fluxi]);

        writeTabbed(os, fluxFieldName);
    }

    os  << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalFlux::sectionalFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sampleFlux(name, runTime, dict)
{
    const fixedSectionalSystem& system =
        mesh_.lookupObject<fixedSectionalSystem>("fixedSectionalSystem");

    fluxFieldNames_.setSize(system.distribution().size());

    forAll(system.distribution(), i)
    {
        const section& sec = system.distribution()[i];

        fluxFieldNames_[i] =
            IOobject::groupName("phiEff", sec.M().name());
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalFlux::~sectionalFlux()
{}

// ************************************************************************* //
