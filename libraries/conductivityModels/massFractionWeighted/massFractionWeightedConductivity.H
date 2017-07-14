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

\file massFractionWeightedConductivity.H
\brief Simple mass-fraction weighted mixture heat conductivity.

This model takes all individual conductivities and computes the mixture
conductivity by a mass fraction weighting.

*/

#ifndef massFractionWeightedConductivity_H
#define massFractionWeightedConductivity_H

#include "conductivityModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace conductivityModels
{

/*---------------------------------------------------------------------------*\
              Class massFractionWeightedConductivity Declaration
\*---------------------------------------------------------------------------*/

class massFractionWeightedConductivity
:
    public conductivityModel
{
    // Private data


    // Private Member Functions

        //- Disallow copy construct
        massFractionWeightedConductivity(const massFractionWeightedConductivity&);

        //- Disallow default bitwise assignment
        void operator=(const massFractionWeightedConductivity&);


protected:


public:

    //- Runtime type information
    TypeName("massFractionWeightedConductivity");


    // Constructors

        //- Construct from components
        massFractionWeightedConductivity
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~massFractionWeightedConductivity();


    // Member Functions

        // Thermodynamics

            //- Update effective conductivity
            void update();

        // I-O

            //- Read info from dictionaries
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace conductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
