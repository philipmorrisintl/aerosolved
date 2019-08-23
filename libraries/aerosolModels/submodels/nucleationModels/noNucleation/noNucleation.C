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

#include "noNucleation.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noNucleation, 0);
addToRunTimeSelectionTable(nucleationModel, noNucleation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noNucleation::noNucleation
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    nucleationModel(type(), aerosol, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noNucleation::~noNucleation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

nucData noNucleation::rate
(
    const scalar& p,
    const scalar& T,
    const scalarList& Y,
    const scalarList& pSat,
    const scalarList& D,
    const scalarList& rhoDisp,
    const scalarList& sigma
) const
{
    return nucData(pSat.size());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
