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

#ifndef aerosolPhase_H
#define aerosolPhase_H

#include "speciesTable.H"
#include "aerosolThermoReader.H"
#include "multiComponentPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class aerosolPhase Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class aerosolPhase
:
    public speciesTable,
    public aerosolThermoReader<ThermoType>,
    public multiComponentPhaseMixture<ThermoType>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        aerosolPhase(const aerosolPhase&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolPhase&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        aerosolPhase(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~aerosolPhase()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return speciesTable::size();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "aerosolPhase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
