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

\file setDeltaT.H
\brief Set the time step size based on three CFL numbers

This file is a copy of the default setDeltaT.H file from OpenFOAM, but also
includes now two additional CFL numbers: a drift CFL and diffusional CFL
number. The time step is limited by the maximum of the three CFL numbers.

*/

if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo/(max(max(max(CoNum, driftCoNum), diffMCoNum), diffCorrYCoNum) + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
