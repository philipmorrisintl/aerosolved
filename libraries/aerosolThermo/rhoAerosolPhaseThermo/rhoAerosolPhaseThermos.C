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

#include "rhoAerosolPhaseThermo.H"
#include "heAerosolRhoThermo.H"
#include "heRhoThermo.H"

#include "makeAerosolThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "thermoPhysicsTypes.H"
#include "dispersedThermoPhysicsTypes.H"
#include "aerosolPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Continuous models

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constGasEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    gasEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleGasEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    incompressibleGasEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    icoPoly8EThermoPhysics
);

// Dispersed models

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constDispEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    dispEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleDispEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressibleFuncDispEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    constIncompressiblePoly8DispEThermoPhysics
);

makeAerosolThermo
(
    rhoThermo,
    rhoAerosolPhaseThermo,
    heAerosolRhoThermo,
    heRhoThermo,
    aerosolPhase,
    adiabaticDispEThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
