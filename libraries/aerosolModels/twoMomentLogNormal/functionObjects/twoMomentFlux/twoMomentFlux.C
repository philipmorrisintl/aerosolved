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

#include "twoMomentFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(twoMomentFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        twoMomentFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::twoMomentFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Two-moment number flux [#/s]");
    writeCommented(os, "Time");

    forAll(fluxFieldNames_, fluxi)
    {
        const word fluxFieldName(fluxFieldNames_[fluxi]);

        writeTabbed(os, fluxFieldName);
    }

    os  << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::twoMomentFlux::twoMomentFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sampleFlux(name, runTime, dict)
{
    fluxFieldNames_.setSize(1);

    fluxFieldNames_[0] = IOobject::groupName("phiEff", "M");

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::twoMomentFlux::~twoMomentFlux()
{}

// ************************************************************************* //
