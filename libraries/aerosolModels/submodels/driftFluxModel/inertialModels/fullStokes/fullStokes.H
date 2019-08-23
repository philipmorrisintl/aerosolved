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

#ifndef fullStokes_H
#define fullStokes_H

#include "inertialModel.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fullStokes Declaration
\*---------------------------------------------------------------------------*/

class fullStokes
:
    public inertialModel
{
private:

    // Private Data

        //- Drift velocity base field
        volVectorField V_;

        //- Relative velocity fields
        HashPtrTable<volVectorField> fields_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- Solver tolerance
        scalar TOL_;


    // Private member functions

        //- Disallow default bitwise copy construct
        fullStokes(const fullStokes&);

        //- Disallow default bitwise assignment
        void operator=(const fullStokes&);

        //- Return the relative velocity field
        volVectorField& Vfield(const word sizeName);


public:

    //- Runtime type information
    TypeName("fullStokes");


    // Constructors

        //- Construct from aerosol model
        fullStokes(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~fullStokes();


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
