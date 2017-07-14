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

/**

\file massFractionWeightedViscosity.H
\brief Simple mass-fraction weighted mixture viscosity.

This model takes all individual viscosities and computes the mixture viscosity
by a mass fraction weighting.

*/

#ifndef massFractionWeightedViscosity_H
#define massFractionWeightedViscosity_H

#include "viscosityModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
              Class massFractionWeightedViscosity Declaration
\*---------------------------------------------------------------------------*/

class massFractionWeightedViscosity
:
    public viscosityModel
{
    // Private data


    // Private Member Functions

        //- Disallow copy construct
        massFractionWeightedViscosity(const massFractionWeightedViscosity&);

        //- Disallow default bitwise assignment
        void operator=(const massFractionWeightedViscosity&);


protected:


public:

    //- Runtime type information
    TypeName("massFractionWeightedViscosity");


    // Constructors

        //- Construct from components
        massFractionWeightedViscosity
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~massFractionWeightedViscosity();


    // Member Functions

        // Thermodynamics

            //- Update effective viscosity
            void update();

        // I-O

            //- Read info from dictionaries
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
