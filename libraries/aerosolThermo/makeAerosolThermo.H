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

#ifndef makeAerosolThermo_H
#define makeAerosolThermo_H

#include "addToRunTimeSelectionTable.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAerosolThermo(BaseThermo,CThermo,AeroPhaseThermo,PhaseThermo,Phase,ThermoPhys)     \
                                                                               \
typedef AeroPhaseThermo                                                        \
<                                                                              \
    CThermo,                                                                   \
    SpecieMixture                                                              \
    <                                                                          \
        Phase                                                                  \
        <                                                                      \
            ThermoPhys                                                         \
        >                                                                      \
    >                                                                          \
> AeroPhaseThermo##Phase##ThermoPhys;                                          \
                                                                               \
typedef PhaseThermo                                                            \
<                                                                              \
    CThermo,                                                                   \
    SpecieMixture                                                              \
    <                                                                          \
        Phase                                                                  \
        <                                                                      \
            ThermoPhys                                                         \
        >                                                                      \
    >                                                                          \
> PhaseThermo##Phase##ThermoPhys;                                              \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    AeroPhaseThermo##Phase##ThermoPhys,                                        \
    (#AeroPhaseThermo"<"#Phase"<" + ThermoPhys::typeName() + ">>").c_str(),    \
    0                                                                          \
);                                                                             \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    PhaseThermo##Phase##ThermoPhys,                                            \
    (#PhaseThermo"<"#Phase"<" + ThermoPhys::typeName() + ">>").c_str(),        \
    0                                                                          \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    AeroPhaseThermo##Phase##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    fluidThermo,                                                               \
    AeroPhaseThermo##Phase##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    AeroPhaseThermo##Phase##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    CThermo,                                                                   \
    AeroPhaseThermo##Phase##ThermoPhys,                                        \
    fvMesh                                                                     \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
