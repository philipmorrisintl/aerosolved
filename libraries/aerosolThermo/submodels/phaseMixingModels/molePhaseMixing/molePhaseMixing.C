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

#include "molePhaseMixing.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(molePhaseMixing, 0);
addToRunTimeSelectionTable(phaseMixingModel, molePhaseMixing, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

molePhaseMixing::molePhaseMixing
(
    const aerosolThermo& thermo
)
:
    phaseMixingModel(type(), thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

molePhaseMixing::~molePhaseMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> molePhaseMixing::mix
(
    const volScalarField& phiCont,
    const volScalarField& phiDisp
) const
{
    return
        thermo_.WMix()
      * (
            (
                thermo_.sumY()
              / thermo_.thermoCont().WMix()
              * phiCont
            )
          + (
                thermo_.sumZ()
              / thermo_.thermoDisp().WMix()
              * phiDisp
            )
        );
}

tmp<scalarField> molePhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const labelList& cells
) const
{
    return
        thermo_.WMix(cells)
      * (
            (
                thermo_.sumY(cells)
              / thermo_.thermoCont().WMix(cells)
              * phiCont
            )
          + (
                thermo_.sumZ(cells)
              / thermo_.thermoDisp().WMix(cells)
              * phiDisp
            )
        );
}

tmp<scalarField> molePhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const label patchi
) const
{
    return
        thermo_.WMix(patchi)
      * (
            (
                thermo_.sumY(patchi)
              / thermo_.thermoCont().WMix(patchi)
              * phiCont
            )
          + (
                thermo_.sumZ(patchi)
              / thermo_.thermoDisp().WMix(patchi)
              * phiDisp
            )
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
