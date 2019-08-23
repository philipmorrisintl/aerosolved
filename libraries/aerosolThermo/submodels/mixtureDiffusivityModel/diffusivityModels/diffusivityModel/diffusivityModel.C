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

#include "diffusivityModel.H"
#include "aerosolThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(diffusivityModel, 0);
    defineRunTimeSelectionTable(diffusivityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::diffusivityModel::diffusivityModel
(
    const word& entryName,
    aerosolThermo& thermo,
    const label j,
    const label k
)
:
    name_(entryName),
    thermo_(thermo),
    j_(j),
    k_(k)
{}

Foam::diffusivityModel::diffusivityModel(const diffusivityModel& de)
:
    refCount(),
    name_(de.name_),
    thermo_(de.thermo_),
    j_(de.j_),
    k_(de.k_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModel::~diffusivityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::diffusivityModel::name() const
{
    return name_;
}

// ************************************************************************* //
