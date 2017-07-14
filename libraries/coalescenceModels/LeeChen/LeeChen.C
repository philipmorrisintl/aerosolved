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

#include "coalescenceModel.H"
#include "LeeChen.H"

#include "makeCoalescenceTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    makeCoalescenceTypes(LeeChen, coalescenceModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::LeeChen::LeeChen
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    coalescenceModel(mesh, aerosol),
    muEffPtr_(NULL)
{
    read();

    const volScalarField& muEff = mesh.lookupObject<volScalarField>("muEff");

    muEffPtr_ = &muEff;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::LeeChen::~LeeChen()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::coalescenceModels::LeeChen::getCoaRateCell
(
    const Foam::scalar vi,
    const Foam::scalar vj,
    const Foam::label jCell,
    const Foam::scalar Kn
)
{
    const volScalarField& T = thermo_.T();
    const volScalarField& muEff = *muEffPtr_;

    scalar K = 2.0 * k() * T[jCell] / (3.0 * muEff[jCell]);

    scalar C = 1.0+1.246*Kn;

    scalar a = pow(vi, 1.0/3.0) + pow(vj, 1.0/3.0);

    return K * a * C * (pow(vi, -1.0/3.0) + pow(vj, -1.0/3.0));
}

bool Foam::coalescenceModels::LeeChen::read()
{
    if (coalescenceModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline const Foam::scalar& Foam::coalescenceModels::LeeChen::W() const
{
    return W_;
}


