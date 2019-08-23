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

#include "dispersedPhase.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dispersedPhase, 0);
addToRunTimeSelectionTable(phaseMixingModel, dispersedPhase, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dispersedPhase::dispersedPhase
(
    const aerosolThermo& thermo
)
:
    phaseMixingModel(type(), thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dispersedPhase::~dispersedPhase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> dispersedPhase::mix
(
    const volScalarField& phiCont,
    const volScalarField& phiDisp
) const
{
    return phiDisp*1.0;
}

tmp<scalarField> dispersedPhase::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const labelList& cells
) const
{
    return phiDisp*1.0;
}

tmp<scalarField> dispersedPhase::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const label patchi
) const
{
    return phiDisp*1.0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
