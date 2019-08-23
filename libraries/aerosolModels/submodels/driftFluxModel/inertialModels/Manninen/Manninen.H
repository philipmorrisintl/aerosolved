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

#ifndef Manninen_H
#define Manninen_H

#include "inertialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Manninen Declaration
\*---------------------------------------------------------------------------*/

class Manninen
:
    public inertialModel
{
private:

    // Private Data

        //- Drift velocity base field
        volVectorField V_;


    // Private member functions

        //- Disallow default bitwise copy construct
        Manninen(const Manninen&);

        //- Disallow default bitwise assignment
        void operator=(const Manninen&);


public:

    //- Runtime type information
    TypeName("Manninen");


    // Constructors

        //- Construct from aerosol model
        Manninen(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~Manninen();


    // Member Functions

        //- Compute the relative velocity given a droplet size field and size
        //  name
        tmp<volVectorField> V(const volScalarField& d, const word sizeName);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
