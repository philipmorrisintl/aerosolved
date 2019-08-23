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

#ifndef aerosolReader_H
#define aerosolReader_H

#include "aerosolReader.H"
#include "fileName.H"
#include "dictionary.H"
#include "speciesTable.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class aerosolReader Declaration
\*---------------------------------------------------------------------------*/

class aerosolReader
{

protected:

    // Protected Member Objects

        //- Phase name
        word phaseName_;

        //- Thermo properties dictionary
        const dictionary& thermoDict_;

        //- Species properties dictionary
        const dictionary& speciesDict_;

        //- Aerosol thermo properties dictionary
        const dictionary aerosolThermoDict_;

        //- Table of species
        speciesTable& speciesTable_;


private:

    // Private Member Functions

        //- Set the species list
        speciesTable& setSpecies
        (
            const dictionary& dict,
            speciesTable& species,
            const word phaseName
        );

        //- Disallow default bitwise copy construct
        aerosolReader(const aerosolReader&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolReader&);


public:

    //- Default aerosol species properties dictionary name
    static const word speciesPropertiesName;

    // Constructors

        //- Construct by getting the foamAerosol and thermodynamics file names
        //  from dictionary
        aerosolReader
        (
            const dictionary& thermoDict,
            speciesTable& species,
            const word& phaseName
        );


    //- Destructor
    virtual ~aerosolReader()
    {}


    // Member functions

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Thermo properties dictionary
        const dictionary& thermoDict() const
        {
            return thermoDict_;
        }

        //- Species properties dictionary
        const dictionary& speciesDict() const
        {
            return speciesDict_;
        }

        //- Aerosol thermo properties dictionary
        const dictionary& aerosolThermoDict() const
        {
            return aerosolThermoDict_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
