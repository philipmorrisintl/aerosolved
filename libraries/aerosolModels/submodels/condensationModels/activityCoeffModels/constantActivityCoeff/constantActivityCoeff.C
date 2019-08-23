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

#include "constantActivityCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantActivityCoeff, 0);
addToRunTimeSelectionTable
(
    activityCoeffModel,
    constantActivityCoeff,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantActivityCoeff::constantActivityCoeff
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    activityCoeffModel(type(), aerosol, dict),
    coeffs_(aerosol.thermo().activeSpecies().size(), 1.0)
{
    forAll(aerosol.thermo().activeSpecies(), j)
    {
        coeffs_[j] = dict.lookupOrDefault<scalar>
            (
                aerosol.thermo().activeSpecies()[j],
                1.0
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantActivityCoeff::~constantActivityCoeff()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalarList constantActivityCoeff::activity(const scalarList& Z) const
{
    return coeffs_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
