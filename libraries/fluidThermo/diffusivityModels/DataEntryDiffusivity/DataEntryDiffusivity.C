/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2017 Philip Morris International

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

#include "DataEntryDiffusivity.H"
#include "makeDiffusivityModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    makeDiffusivityModels(DataEntryDiffusivity, diffusivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::DataEntryDiffusivity::DataEntryDiffusivity
(
    const word& entryName,
    const dictionary& dict,
    const dictionary species,
    const label a,
    const label b
)
:
    diffusivityModel(entryName, species, a, b),
    diffusivity_(DataEntry<scalar>::New(entryName, dict))
{
}


Foam::diffusivityModels::DataEntryDiffusivity::DataEntryDiffusivity
(
    const DataEntryDiffusivity& fsg
)
:
    diffusivityModel(fsg),
    diffusivity_(fsg.diffusivity_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModels::DataEntryDiffusivity::~DataEntryDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::diffusivityModels::DataEntryDiffusivity::value
(
    const Foam::scalar T,
    const Foam::scalar p
) const
{
    return diffusivity_().value(T);
}

// ************************************************************************* //
