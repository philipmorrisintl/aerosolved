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

#ifndef noInertial_H
#define noInertial_H

#include "inertialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noInertial Declaration
\*---------------------------------------------------------------------------*/

class noInertial
:
    public inertialModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noInertial(const noInertial&);

        //- Disallow default bitwise assignment
        void operator=(const noInertial&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from aerosol model
        noInertial(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~noInertial();


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
