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

#include "heAerosolRhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heAerosolRhoThermo<BasicPsiThermo, MixtureType>::heAerosolRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heRhoThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{}

template<class BasicPsiThermo, class MixtureType>
Foam::heAerosolRhoThermo<BasicPsiThermo, MixtureType>::heAerosolRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    heRhoThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictName)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heAerosolRhoThermo<BasicPsiThermo, MixtureType>::~heAerosolRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heAerosolRhoThermo<BasicPsiThermo, MixtureType>::WMix
(
    const Foam::labelList& cells
) const
{
    tmp<scalarField> tW(new scalarField(cells.size(), 0.0));

    scalarField& W = tW.ref();

    forAll(W, i)
    {
        W[i] = this->cellMixture(cells[i]).W();
    }

    return tW;
}

template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heAerosolRhoThermo<BasicPsiThermo, MixtureType>::WMix
(
    const Foam::label patchi
) const
{
    const fvMesh& mesh = this->T_.mesh();

    const polyPatch& patch = mesh.boundaryMesh()[patchi];

    tmp<scalarField> tW(new scalarField(patch.size(), 0.0));

    scalarField& W = tW.ref();

    forAll(patch, facei)
    {
        W[facei] = this->patchFaceMixture(patchi, facei).W();
    }

    return tW;
}

// ************************************************************************* //
