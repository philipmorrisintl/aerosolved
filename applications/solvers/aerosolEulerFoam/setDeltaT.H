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

if (adjustTimeStep)
{
    scalar maxConvDeltaTFact = maxCo/(CoNum + SMALL);
    scalar maxAerosolDeltaTFact = maxAerosolCo/(aerosolCoNum + SMALL);
    scalar maxDriftDeltaTFact = maxDriftCo/(driftCoNum + SMALL);

    scalar maxDeltaTFact
    (
        min
        (
            min(maxConvDeltaTFact, maxAerosolDeltaTFact),
            maxDriftDeltaTFact
        )
    );

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
