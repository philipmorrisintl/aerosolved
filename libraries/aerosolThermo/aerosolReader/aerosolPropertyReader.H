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

#ifndef aerosolPropertyReader_H
#define aerosolPropertyReader_H

#include "aerosolReader.H"
#include "Function1.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class aerosolPropertyReader Declaration
\*---------------------------------------------------------------------------*/

class aerosolPropertyReader
:
    public aerosolReader
{

        HashPtrTable<Function1<scalar>, word> properties_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        aerosolPropertyReader(const aerosolPropertyReader&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolPropertyReader&);


public:

    // Constructors

        //- Construct by getting the foamAerosol and thermodynamics file names
        //  from dictionary
        aerosolPropertyReader
        (
            const dictionary& thermoDict,
            speciesTable& species,
            const word& phaseName
        );

        //- Check if the property for the given species exists
        bool propertyFound
        (
            const word& speciesName,
            const word& propertyName
        ) const;

        //- Return the property function by name
        const Function1<scalar>& property
        (
            const word& speciesName,
            const word& propertyName
        );

        //- Return the property function by index
        const Function1<scalar>& property
        (
            const label& j,
            const word& propertyName
        );


    //- Destructor
    virtual ~aerosolPropertyReader()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
