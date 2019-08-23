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

#include "volumePhaseMixing.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volumePhaseMixing, 0);
addToRunTimeSelectionTable(phaseMixingModel, volumePhaseMixing, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volumePhaseMixing::volumePhaseMixing
(
    const aerosolThermo& thermo
)
:
    phaseMixingModel(type(), thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volumePhaseMixing::~volumePhaseMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> volumePhaseMixing::mix
(
    const volScalarField& phiCont,
    const volScalarField& phiDisp
) const
{
    return
        thermo_.rho()
      * (
            (
                thermo_.sumY()
              / thermo_.thermoCont().rho()
              * phiCont
            )
          + (
                thermo_.sumZ()
              / thermo_.thermoDisp().rho()
              * phiDisp
            )
        );
}

tmp<scalarField> volumePhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const labelList& cells
) const
{
    return
        thermo_.rho(cells)
      * (
            (
                thermo_.sumY(cells)
              / thermo_.thermoCont().rho(cells)
              * phiCont
            )
          + (
                thermo_.sumZ(cells)
              / thermo_.thermoDisp().rho(cells)
              * phiDisp
            )
        );
}

tmp<scalarField> volumePhaseMixing::mix
(
    const scalarField& phiCont,
    const scalarField& phiDisp,
    const label patchi
) const
{
    return
        thermo_.rho(patchi)
      * (
            (
                thermo_.sumY(patchi)
              / thermo_.thermoCont().rho(patchi)
              * phiCont
            )
          + (
                thermo_.sumZ(patchi)
              / thermo_.thermoDisp().rho(patchi)
              * phiDisp
            )
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
