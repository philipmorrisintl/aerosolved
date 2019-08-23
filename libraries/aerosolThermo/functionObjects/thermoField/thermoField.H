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

/**

\file thermoField.H
\brief functionObject class to write a thermo field to disk

This functionObject class can be used to write additional thermo fields to
disk, such as Cp, Cv or nu.

*/

#ifndef functionObjects_thermoField_H
#define functionObjects_thermoField_H

#include "fvMeshFunctionObject.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class thermoField Declaration
\*---------------------------------------------------------------------------*/

class thermoField
:
    public fvMeshFunctionObject
{
public:

    // Public data types

        //- Thermo field type enumeration
        enum thermoFieldType
        {
            Cp,
            Cv,
            Cpv,
            CpByCpv,
            gamma,
            hc,
            nu,
            kappa,
            sumY,
            sumZ,
            lambda,
        };

        //- Thermo field type names
        static const Enum<thermoFieldType> thermoFieldTypeNames_;


private:

    // Private Data

        //- Thermo field to obtain
        thermoFieldType thermoField_;

        //- Const reference to the aerosol model
        aerosolThermo& thermo_;

        //- Name of the result field
        word resultName_;


    // Private Member Functions

        //- Calculate the diameter field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("thermoField");


    // Constructors

        //- Construct from Time and dictionary
        thermoField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermoField();

    //- Member Functions
        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
