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

#include "Manninen.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Manninen, 0);
addToRunTimeSelectionTable(inertialModel, Manninen, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Manninen::Manninen
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    inertialModel(type(), aerosol, dict),
    V_
    (
        IOobject
        (
            "V",
            aerosol.mesh().time().timeName(),
            aerosol.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aerosol.mesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Manninen::~Manninen()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volVectorField> Manninen::V
(
    const volScalarField& d,
    const word sizeName
)
{
    const fvMesh& mesh = aerosol_.mesh();

    const volScalarField& rho = aerosol_.rho();
    const volScalarField& rhol = aerosol_.thermo().thermoDisp().rho();
    const volScalarField& rhog = aerosol_.thermo().thermoCont().rho();

    const volVectorField& U = aerosol_.U();
    const surfaceScalarField& phi = aerosol_.phi();

    const volScalarField& mug = aerosol_.thermo().thermoCont().mu();

    const volScalarField tau(Foam::sqr(d)*rhol/(18.0*mug));
    const volVectorField G((rhol-rhog)/rhol*g_);

    const surfaceScalarField phiU(phi/linearInterpolate(rho));

    tmp< fv::gaussConvectionScheme<vector> > tdivScheme
    (
        new fv::gaussConvectionScheme<vector>
        (
            mesh,
            phiU,
            upwind<vector>(mesh, phiU)
        )
    );

    fv::gaussConvectionScheme<vector>& divScheme = tdivScheme.ref();

    const volVectorField a
    (
        fvc::ddt(U)
      + divScheme.fvcDiv(phiU, U)
      - fvc::Sp(fvc::div(phiU), U)
    );

    tmp<volVectorField> tV
    (
        new volVectorField
        (
            IOobject
            (
                "V",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (G-a)*tau,
            V_.boundaryField().types()
        )
    );

    volVectorField& V = tV.ref();

    V.correctBoundaryConditions();

    return limit(V);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
