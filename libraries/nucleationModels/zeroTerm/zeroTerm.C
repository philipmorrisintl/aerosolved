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

#include "nucleationModel.H"
#include "zeroTerm.H"

#include "makeNucleationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{
    makeNucleationTypes(zeroTerm, nucleationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    nucleationModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField>&
Foam::nucleationModels::zeroTerm::getNucFields()
{
    label n = thermo_.nSpeciesPhaseChange();

    forAll(thermo_.speciesPhaseChange(), j)
    {
        SJDnuc_[j].internalField() = 0.0;
    }

    SJDnuc_[n].internalField() = 0.0;
    SJDnuc_[n+1].internalField() = 0.0;

    return SJDnuc_;
}

bool Foam::nucleationModels::zeroTerm::read()
{
    if (nucleationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
