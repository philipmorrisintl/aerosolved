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

#ifndef aerosolSubModelBase_H
#define aerosolSubModelBase_H

#include "subModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward decleration
class aerosolModel;

/*---------------------------------------------------------------------------*\
                      Class aerosolSubModelBase Declaration
\*---------------------------------------------------------------------------*/

class aerosolSubModelBase
:
    public subModelBase
{
protected:

    // Protected data

        //- Reference to aerosol model
        aerosolModel& aerosol_;


public:

    // Constructors

        //- Construct null
        aerosolSubModelBase(aerosolModel& aerosol);

        //- Construct from aerosol model
        aerosolSubModelBase
        (
            aerosolModel& aerosol,
            const dictionary& dict,
            const word& baseName,
            const word& modelType,
            const word& dictExt = "Coeffs"
        );

        //- Construct from aerosol model
        aerosolSubModelBase
        (
            const word& modelName,
            aerosolModel& aerosol,
            const dictionary& dict,
            const word& baseName,
            const word& modelType
        );


    //- Destructor
    virtual ~aerosolSubModelBase();


    // Member Functions

        // Access

            //- Flag to indicate when to write a property
            virtual bool writeTime() const;

            //- Return const access to the aerosol model
            inline const aerosolModel& aerosol() const;

            //- Return the reference to the aerosol model
            inline aerosolModel& aerosol();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aerosolSubModelBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
