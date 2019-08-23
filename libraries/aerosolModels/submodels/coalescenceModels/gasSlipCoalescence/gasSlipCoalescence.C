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

#include "gasSlipCoalescence.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gasSlipCoalescence, 0);
addToRunTimeSelectionTable(coalescenceModel, gasSlipCoalescence, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gasSlipCoalescence::gasSlipCoalescence
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    coalescenceModel(type(), aerosol, dict),
    A_(readScalar(dict.lookup("A"))),
    p_(4),
    q_(4)
{
    p_[0] = 0.0;
    p_[1] = 1.0;
    p_[2] = 0.0;
    p_[3] = 1.0;

    q_[0] = 0.0;
    q_[1] = -1.0;
    q_[2] = -1.0;
    q_[3] = -2.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gasSlipCoalescence::~gasSlipCoalescence()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

coaData gasSlipCoalescence::rate
(
    const scalar& p,
    const scalar& T,
    const scalar& mu,
    const scalar& rhog,
    const scalar& rhol,
    const scalar& d
) const
{
    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar R = constant::physicoChemical::R.value();

    const scalar K(2.0*kB*T/(3.0*mu));

    const scalar m(rhog*R*T/p/NA);

    const scalar lambda(mu/p*sqrt(pi*kB*T/(2.0*m)));

    scalarList w(4, K);

    w[2] *= A_*lambda*2.0;
    w[3] *= A_*lambda*2.0;

    return coaData(w, p_, q_, true);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
