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

#include "phaseMixing.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseMixing::phaseMixing
(
    const aerosolThermo& thermo
)
:
    thermo_(thermo),
    viscosity_(),
    conductivity_(),
    heatCapacity_()
{
    viscosity_ =
        phaseMixingModel::New
        (
            thermo,
            word(thermo.subDict("phaseMixing").lookup("viscosity"))
        );

    conductivity_ =
        phaseMixingModel::New
        (
            thermo,
            word(thermo.subDict("phaseMixing").lookup("conductivity"))
        );

    heatCapacity_ =
        phaseMixingModel::New
        (
            thermo,
            word(thermo.subDict("phaseMixing").lookup("heatCapacity"))
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

phaseMixing::~phaseMixing()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
