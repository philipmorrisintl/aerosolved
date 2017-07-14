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

#include "zeroTerm.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
    makeFluidThermoTypes(zeroTerm, aerosolModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::zeroTerm::zeroTerm
(
    const fvMesh& mesh
)
:
    aerosolModel(mesh),
    dropletSizeDimension_(dimless),
    mesh_(mesh)
{
    this->read();
    updateSizeDistribution();

    if(sizeDistType_ != NOSIZEDIST)
    {
        FatalErrorIn("Foam::aerosolModels::zeroTerm::zeroTerm(const fvMesh& mesh)")
            << "This is a moment method. The size distribution type must be 'none'." << exit(FatalError);
    }

    M_.setSize(1);
    J_.setSize(1);
    S_.setSize(thermo().nSpecies());
    phid_.setSize(1);

    M_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "M",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    J_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "J",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("J", dimless/dimVolume/dimTime, 0.0)
        )
    );

    forAll(thermo().species(), j)
    {
        S_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word("S." + Foam::name(j)),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("S", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }

    phid_.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                "phid",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("phid", dimVelocity*dimDensity*dimArea, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::zeroTerm::~zeroTerm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::zeroTerm::update()
{
}

void Foam::aerosolModels::zeroTerm::fractionalStepInternal()
{
}

void Foam::aerosolModels::zeroTerm::fractionalStepExternal()
{
}

void Foam::aerosolModels::zeroTerm::checkConsistency()
{
}

void Foam::aerosolModels::zeroTerm::correctSizeDistribution()
{
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::zeroTerm::dcm()
{
    tmp<volScalarField> tdcm
    (
        new volScalarField
        (
            IOobject
            (
                "dcm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dcm", dimLength, 0.0)
        )
    );

    //FatalErrorIn("Foam::aerosolModels::zeroTerm::dcm()")
    //    << "Not yet implemented." << exit(FatalError);

    return tdcm;
}

Foam::tmp<Foam::volScalarField> Foam::aerosolModels::zeroTerm::dmm()
{
    tmp<volScalarField> tdmm
    (
        new volScalarField
        (
            IOobject
            (
                "dmm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dmm", dimLength, 0.0)
        )
    );

    //FatalErrorIn("Foam::aerosolModels::zeroTerm::dmm()")
    //    << "Not yet implemented." << exit(FatalError);

    return tdmm;
}

bool Foam::aerosolModels::zeroTerm::read()
{
    if (aerosolModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
