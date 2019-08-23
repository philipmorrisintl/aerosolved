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

#include "FullerSchettlerGiddings.H"
#include "aerosolThermo.H"
#include "makeDiffusivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    makeDiffusivityModel(FullerSchettlerGiddings, diffusivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::FullerSchettlerGiddings::FullerSchettlerGiddings
(
    const word& entryName,
    const dictionary& dict,
    aerosolThermo& thermo,
    const label j,
    const label k
)
:
    diffusivityModel(entryName, thermo, j, k),
    Vdj_(thermo.thermoCont().property(j, "Vd")),
    Vdk_(thermo.thermoCont().property(k, "Vd"))
{}


Foam::diffusivityModels::FullerSchettlerGiddings::FullerSchettlerGiddings
(
    const FullerSchettlerGiddings& se
)
:
    diffusivityModel(se),
    Vdj_(se.Vdj_),
    Vdk_(se.Vdk_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModels::FullerSchettlerGiddings::~FullerSchettlerGiddings()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::diffusivityModels::FullerSchettlerGiddings::D() const
{
    const scalarField& T = thermo_.thermoCont().T().field();
    const scalarField& p = thermo_.thermoCont().p().field();

    const scalar Wj(thermo_.thermoCont().composition().W(j_));
    const scalar Wk(thermo_.thermoCont().composition().W(k_));

    return
        1E-7 * pow(T, 1.75) / (p/1.01325E5) * Foam::sqrt(1.0/Wj + 1.0/Wk)
      / sqr(pow(Vdj_.value(T), 1.0/3.0) + pow(Vdk_.value(T), 1.0/3.0));
}

// ************************************************************************* //
