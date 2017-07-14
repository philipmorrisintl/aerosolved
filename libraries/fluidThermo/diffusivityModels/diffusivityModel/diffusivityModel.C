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

#include "diffusivityModel.H"

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
    const dictionary species,
    const label a,
    const label b
)
:
    refCount(),
    name_(entryName),
    species_(species),
    a_(a),
    b_(b)
{}


Foam::diffusivityModel::diffusivityModel(const diffusivityModel& de)
:
    refCount(),
    name_(de.name_),
    species_(de.species_),
    a_(de.a_),
    b_(de.b_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModel::~diffusivityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::diffusivityModel::name() const
{
    return name_;
}

namespace Foam // Necessary for Doxygen
{

scalar diffusivityModel::value
(
    const scalar T,
    const scalar p
) const
{
    notImplemented("Type Foam::diffusivityModel::value(const scalar, const scalar) const");

    return scalar(0.0);
}


tmp<scalarField> diffusivityModel::value
(
    const scalarField& T,
    const scalarField& p
) const
{
    tmp<scalarField> tfld(new scalarField(T.size()));
    scalarField& fld = tfld();

    forAll(T, i)
    {
        fld[i] = this->value(T[i], p[i]);
    }
    return tfld;
}

} // End Foam namespace

// ************************************************************************* //
