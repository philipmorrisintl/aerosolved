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

#include "function.H"
#include "aerosolThermo.H"
#include "makeDiffusivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    makeDiffusivityModel(function, diffusivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::diffusivityModels::function::read(const dictionary& coeffs)
{}


Foam::diffusivityModels::function::function
(
    const word& entryName,
    const dictionary& dict,
    aerosolThermo& thermo,
    const label j,
    const label k
)
:
    diffusivityModel(entryName, thermo, j, k),
    function_(Function1<scalar>::New(entryName, dict))
{
    read(dict);
}


Foam::diffusivityModels::function::function
(
    const function& se
)
:
    diffusivityModel(se),
    function_(se.function_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModels::function::~function()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::diffusivityModels::function::D() const
{
    return function_->value(thermo_.thermoCont().T().field());
}

// ************************************************************************* //
