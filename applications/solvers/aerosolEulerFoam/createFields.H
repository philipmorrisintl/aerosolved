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

#include "createDriftRDeltaT.H"

Info<< "Creating aerosol model\n" << endl;

autoPtr<aerosolModel> aerosol(aerosolModel::New(mesh));

aerosolThermo& thermo = aerosol->thermo();
fluidThermo& fthermo = aerosol->thermo();

rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();
rhoAerosolPhaseThermo& thermoDisp = thermo.thermoDisp();

basicMultiComponentMixture& compCont = thermoCont.composition();
basicMultiComponentMixture& compDisp = thermoDisp.composition();

PtrList<volScalarField>& Y = compCont.Y();
PtrList<volScalarField>& Z = compDisp.Y();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info<< "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        fthermo
    )
);

// Set the turbulence into the aerosol model
aerosol->setTurbulence(turbulence());


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}

forAll(Z, i)
{
    fields.add(Z[i]);
}

#include "createMRF.H"

PtrList<surfaceScalarField> phiEffY(Y.size());
PtrList<surfaceScalarField> phiEffZ(Z.size());

forAll(Y, j)
{
    mesh.setFluxRequired(Y[j].name());

    phiEffY.set
    (
        j,
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("phiEff", Y[j].name()),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi*0.0
        )
    );
}

forAll(Z, j)
{
    mesh.setFluxRequired(Z[j].name());

    phiEffZ.set
    (
        j,
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("phiEff", Z[j].name()),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi*0.0
        )
    );
}

scalar calcAerosolAfter=runTime.startTime().value();

if(LTS && !pimple.frozenFlow())
{
    calcAerosolAfter=readScalar(pimple.dict().lookup("calcAerosolAfter"));

    Info<< "Starting calculation of the aerosol after t = "
        << calcAerosolAfter << endl;
}
