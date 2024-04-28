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

#include "rhoAerosolPhaseThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoAerosolPhaseThermo, 0);
    defineRunTimeSelectionTable(rhoAerosolPhaseThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoAerosolPhaseThermo::rhoAerosolPhaseThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    speciesTable(),
    rhoReactionThermo(mesh, phaseName),
    aerosolPropertyReader(*this, *this, phaseName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoAerosolPhaseThermo> Foam::rhoAerosolPhaseThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoAerosolPhaseThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rhoAerosolPhaseThermo::WMix() const
{
    return W();
}


Foam::tmp<Foam::scalarField> Foam::rhoAerosolPhaseThermo::rho
(
    const labelList& cells
) const
{
    tmp<scalarField> trho(new scalarField(cells.size(), 0.0));

    scalarField& rho = trho.ref();

    forAll(cells, i)
    {
        rho[i] = rho_[cells[i]];
    }

    return trho;
}


// ************************************************************************* //
