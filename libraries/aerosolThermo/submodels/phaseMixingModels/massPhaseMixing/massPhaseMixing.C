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

#include "massPhaseMixing.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(massPhaseMixing, 0);
addToRunTimeSelectionTable(phaseMixingModel, massPhaseMixing, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

massPhaseMixing::massPhaseMixing
(
    const aerosolThermo& thermo
)
:
    phaseMixingModel(type(), thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massPhaseMixing::~massPhaseMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> massPhaseMixing::mix
(
    const volScalarField& phiCont,
    const volScalarField& phiDisp
) const
{
    return thermo_.sumY()*phiCont + thermo_.sumZ()*phiDisp;
}

tmp<scalarField> massPhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const labelList& cells
) const
{
    return
        thermo_.sumY(cells)*phiCont
      + thermo_.sumZ(cells)*phiDisp;
}

tmp<scalarField> massPhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const label patchi
) const
{
    return
        thermo_.sumY(patchi)*phiCont
      + thermo_.sumZ(patchi)*phiDisp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
