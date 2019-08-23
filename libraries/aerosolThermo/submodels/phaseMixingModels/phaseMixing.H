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

\file phaseMixing.H
\brief Class which contains the three phaseMixingModels

The phaseMixingModels are for viscosity, conductivity and heat capacity

*/

#ifndef phaseMixing_H
#define phaseMixing_H

#include "autoPtr.H"
#include "phaseMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class phaseMixing Declaration
\*---------------------------------------------------------------------------*/

class phaseMixing
{
private:

    // Private Data

        //- Reference to the aerosol thermo object
        const aerosolThermo& thermo_;

        //- Viscosity mixing object
        autoPtr<phaseMixingModel> viscosity_;

        //- Conductivity mixing object
        autoPtr<phaseMixingModel> conductivity_;

        //- Heat capacity mixing object
        autoPtr<phaseMixingModel> heatCapacity_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseMixing(const phaseMixing&);

        //- Disallow default bitwise assignment
        void operator=(const phaseMixing&);


public:

    // Constructors

        //- Construct from type name, dictionary and thermo model
        phaseMixing(const aerosolThermo& thermo);


    //- Destructor
    virtual ~phaseMixing();


    // Member Functions

        // Access

            //- Reference to the aerosol thermo object
            inline const aerosolThermo& thermo()
            {
                return thermo_;
            }

        // Mixing

            //- Mixing of viscosity
            inline const phaseMixingModel& viscosity() const
            {
                return viscosity_();
            }

            //- Mixing of conductivity
            inline const phaseMixingModel& conductivity() const
            {
                return conductivity_();
            }

            //- Mixing of heat capacity
            inline const phaseMixingModel& heatCapacity() const
            {
                return heatCapacity_();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
