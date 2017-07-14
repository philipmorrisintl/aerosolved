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

#include "condensationEvaporationModel.H"
#include "zeroTerm.H"

#include "makeCondensationEvaporationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{
    makeCondensationEvaporationTypes(zeroTerm, condensationEvaporationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    condensationEvaporationModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::zeroTerm::getCondRateListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{
    const label& nSpeciesPhaseChange = thermo_.nSpeciesPhaseChange();

    List<List<scalar> > I(z.size());

    forAll(z, i)
    {
        I[i] = List<scalar>(nSpeciesPhaseChange, 0.0);
    }

    return I;
}

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::zeroTerm::getEtaGammaListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{
    const label& nSpeciesPhaseChange = thermo_.nSpeciesPhaseChange();

    List< List<scalar> > I(z.size());

    forAll(z, i)
    {
        I[i] = List<scalar>(nSpeciesPhaseChange+1, 0.0);
    }

    return I;
}

Foam::List<Foam::scalar> Foam::condensationEvaporationModels::zeroTerm::psiInv
(
    const Foam::List<Foam::scalar>& zeta,
    const Foam::label jCell
)
{
    List<scalar> a(zeta.size(), -1.0);

    return a;
}

bool Foam::condensationEvaporationModels::zeroTerm::read()
{
    if (condensationEvaporationModel::read())
    {
        // Read condensationEvaporation model coefficients

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
