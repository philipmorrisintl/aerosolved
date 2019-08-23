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

#include "aerosolPropertyReader.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::aerosolPropertyReader::aerosolPropertyReader
(
    const dictionary& thermoDict,
    speciesTable& species,
    const word& phaseName
)
:
    aerosolReader(thermoDict, species, phaseName),
    properties_()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::aerosolPropertyReader::propertyFound
(
    const Foam::word& speciesName,
    const Foam::word& propertyName
) const
{
    const dictionary& properties =
        speciesDict_.subDict(speciesName).subDict("properties");

    return properties.found(propertyName);
}

const Foam::Function1<Foam::scalar>& Foam::aerosolPropertyReader::property
(
    const Foam::word& speciesName,
    const Foam::word& propertyName
)
{
    if (!properties_.found(speciesName+propertyName))
    {
        const dictionary& properties =
            speciesDict_.subDict(speciesName).subDict("properties");

        if (!properties.found(propertyName))
        {
            FatalErrorInFunction
                << "Could not find property " << propertyName
                << " for species " << speciesName << nl
                << abort(FatalError);
        }

        autoPtr<Function1<scalar>> functionPtr =
            Function1<scalar>::New(propertyName, properties);

        properties_.insert(speciesName+propertyName, functionPtr);
    }

    return *properties_[speciesName+propertyName];
}

const Foam::Function1<Foam::scalar>& Foam::aerosolPropertyReader::property
(
    const Foam::label& j,
    const Foam::word& propertyName
)
{
    return property(speciesTable_[j], propertyName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
