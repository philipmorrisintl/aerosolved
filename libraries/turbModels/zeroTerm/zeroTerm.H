/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2017 Philip Morris International

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

#ifndef zeroTerm_H
#define zeroTerm_H

#include "turbModel.H"
#include "DataEntry.H"
#include "fvcGrad.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbModels
{

/*---------------------------------------------------------------------------*\
                            Class zeroTerm Declaration
\*---------------------------------------------------------------------------*/

class zeroTerm
:
    public turbModel
{
    // Private data

    // Private Member Functions

        //- Disallow copy construct
        zeroTerm(const zeroTerm&);

        //- Disallow default bitwise assignment
        void operator=(const zeroTerm&);


protected:


public:

    //- Runtime type information
    TypeName("zeroTerm");


    // Constructors

        //- Construct from components
        zeroTerm
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~zeroTerm();


    // Member Functions

        // Thermodynamics

            //- Update effective turbulence
            void update();

        // I-O

            //- Read info from dictionaries
            bool read();

           tmp<volScalarField> K(volVectorField& U);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
