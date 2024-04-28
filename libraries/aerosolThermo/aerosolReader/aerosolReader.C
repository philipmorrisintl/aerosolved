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

#include "aerosolReader.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::aerosolReader::speciesPropertiesName
(
    "speciesProperties"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::speciesTable& Foam::aerosolReader::setSpecies
(
    const dictionary& dict,
    speciesTable& species,
    const word phaseName
)
{
    wordList sActive(dict.lookup("activeSpecies"));
    wordList sInactive(dict.lookup("inactiveSpecies"));

    word inertSpecie(dict.lookup("inertSpecie"));

    // if (!sInactive.found(inertSpecie))
    if (findIndex(sInactive, inertSpecie) == -1)
    {
        FatalErrorInFunction
            << "The inert specie was not found in the inactive species lists"
            << nl
            << abort(FatalError);
    }

    // if (sActive.found(inertSpecie))
    if (findIndex(sActive, inertSpecie) != -1)
    {
        FatalErrorInFunction
            << "The inert specie cannot be set as active" << nl
            << abort(FatalError);
    }

    forAll(sActive, j)
    {
        // if (sActive.found(sActive[j]))
        if (findIndex(sInactive, sActive[j]) != -1)
        {
            FatalErrorInFunction
                << "Found specie which is both active and inactive ("
                << sActive[j] << ")" << nl
                << abort(FatalError);
        }
    }

    if (phaseName == "dispersed")
    {
        // Order of dispersed species:
        //  1) active species
        //  2) non-inert inactive species

        wordList sl(sActive);

        forAll(sInactive, j)
        {
            if (sInactive[j] != inertSpecie)
            {
                sl.append(sInactive[j]);
            }
        }

        species.transfer(sl);
    }
    else
    {
        // Order of continuous species:
        //  1) active species
        //  2) non-inert inactive species
        //  3) inert inactive specie

        wordList sg(sActive);

        forAll(sInactive, j)
        {
            if (sInactive[j] != inertSpecie)
            {
                sg.append(sInactive[j]);
            }
        }

        sg.append(inertSpecie);

        species.transfer(sg);
    }

    return species;
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::aerosolReader::aerosolReader
(
    const dictionary& thermoDict,
    speciesTable& species,
    const word& phaseName
)
:
    phaseName_(phaseName),
    thermoDict_(thermoDict),
    speciesDict_(thermoDict.subDict("species")),
    aerosolThermoDict_
    (
        IFstream(fileName("constant/thermophysicalProperties").expand())()
    ),
    speciesTable_(setSpecies(aerosolThermoDict_, species, phaseName))
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
