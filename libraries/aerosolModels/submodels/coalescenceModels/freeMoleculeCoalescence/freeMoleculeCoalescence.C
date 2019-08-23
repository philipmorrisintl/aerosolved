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

#include "freeMoleculeCoalescence.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeMoleculeCoalescence, 0);
addToRunTimeSelectionTable
(
    coalescenceModel,
    freeMoleculeCoalescence,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeMoleculeCoalescence::freeMoleculeCoalescence
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    coalescenceModel(type(), aerosol, dict),
    b_(readScalar(dict.lookup("b"))),
    p_(3),
    q_(3)
{
    p_[0] = 0.5;
    p_[1] = 2.0;
    p_[2] = 1.0;

    q_[0] = 0.0;
    q_[1] = -3.0/2.0;
    q_[2] = -0.5;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

freeMoleculeCoalescence::~freeMoleculeCoalescence()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

coaData freeMoleculeCoalescence::rate
(
    const scalar& p,
    const scalar& T,
    const scalar& mu,
    const scalar& rhog,
    const scalar& rhol,
    const scalar& d
) const
{
    const scalar kB = constant::physicoChemical::k.value();

    const scalar K(2.0*kB*T/(3.0*mu));

    // Otto et al. (1994), Eq. (5)

    const scalar Kt(3.0*sqrt(3.0)/2.0*sqrt(sqr(mu)/(rhol*kB*T))*K);

    scalarList w(3, b_*Kt);

    w[2] *= 2.0;

    return coaData(w, p_, q_, true);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
