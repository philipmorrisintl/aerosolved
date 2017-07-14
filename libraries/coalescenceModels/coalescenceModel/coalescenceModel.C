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

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(coalescenceModel, 0);
    defineRunTimeSelectionTable(coalescenceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModel::coalescenceModel
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    IOdictionary
    (
        IOobject
        (
            "coalescenceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    aerosol_(aerosol),
    thermo_(aerosol.thermo()),
    coeffs_(subDict("Coeffs"))
{
    aerosol.getCoaRateField = boost::bind
    (
        &coalescenceModel::getCoaRateField,
        this,
        _1,
        _2,
        _3
    );

    aerosol.getCoaRateCell = boost::bind
    (
        &coalescenceModel::getCoaRateCell,
        this,
        _1,
        _2,
        _3,
        _4
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModel::~coalescenceModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::coalescenceModel::getCoaRateField
(
    const Foam::scalar vi,
    const Foam::scalar vj,
    const Foam::volScalarField& Kn
)
{
    tmp<volScalarField> tBeta
    (
        new volScalarField
        (
            IOobject
            (
                "beta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("beta", dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& beta = tBeta();

    forAll(mesh_.C(), jCell)
    {
        beta[jCell] = getCoaRateCell(vi, vj, jCell, Kn[jCell]);
    }

    return tBeta;
}

bool Foam::coalescenceModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
