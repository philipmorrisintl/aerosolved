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

#ifndef aerosolThermoReader_H
#define aerosolThermoReader_H

#include "HashPtrTable.H"
#include "aerosolReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class aerosolThermoReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class aerosolThermoReader
:
    public aerosolReader
{
        //- Table of the thermodynamic data given in the foamAerosol file
        HashPtrTable<ThermoType> speciesThermo_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        aerosolThermoReader(const aerosolThermoReader&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolThermoReader&);


public:

    // Constructors

        //- Construct by getting the foamAerosol and thermodynamics file names
        //  from dictionary
        aerosolThermoReader
        (
            const dictionary& thermoDict,
            speciesTable& species,
            const word& phaseName
        );


    //- Destructor
    virtual ~aerosolThermoReader()
    {}


    // Member functions

        //- Table of the thermodynamic data given in the foamAerosol file
        const HashPtrTable<ThermoType>& speciesThermo() const
        {
            return speciesThermo_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "aerosolThermoReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
