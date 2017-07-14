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

/**

\file createFields.H
\brief Create fields, model objects and references

In this file a lot of required fields and objects are created, e.g., the
velocity field U, pressure field p and face flux field phi. Also all models
are created in terms of corresponding objects and necessary references
are made to fields contained in these objects.

*/

// Create the thermo object

autoPtr<aerosolModel> tAerosolModel(aerosolModel::New(mesh));
aerosolModel& aerosol = tAerosolModel();

fluidThermo& thermo = aerosol.thermo();

// Create references to thermo object members

PtrList<volScalarField>& Y = thermo.Y();
PtrList<volScalarField>& Z = thermo.Z();

dimensionedScalar p0 = thermo.p0();
volScalarField& p1 = thermo.p1();
volScalarField& psi = thermo.psi();
volScalarField& T = thermo.T();

// List of M fields (only one for moment method, multiple for sectional)

PtrList<volScalarField>& M = aerosol.M();

volScalarField& rho = thermo.rho();
const volScalarField& CpEff = thermo.CpEff();

const dictionary& species = thermo.species();

// Create derived run time selectable models

autoPtr<diffusionModel> tDiffusion(diffusionModel::New(mesh, thermo));
autoPtr<brownianDiffusionModel> tBrownianDiffusion(brownianDiffusionModel::New(mesh, aerosol));
autoPtr<conductivityModel> tConductivity(conductivityModel::New(mesh, thermo));
autoPtr<viscosityModel> tViscosity(viscosityModel::New(mesh, thermo));
autoPtr<nucleationModel> tNucleation(nucleationModel::New(mesh, aerosol));
autoPtr<condensationEvaporationModel> tCondensationEvaporation(condensationEvaporationModel::New(mesh, aerosol));
autoPtr<coalescenceModel> tCoalescence(coalescenceModel::New(mesh, aerosol));
autoPtr<driftVelocityModel> tDriftVelocity(driftVelocityModel::New(mesh, aerosol));


diffusionModel& diffusion = tDiffusion();
brownianDiffusionModel& brownianDiffusion = tBrownianDiffusion();
conductivityModel& conductivity = tConductivity();
viscosityModel& viscosity = tViscosity();

// Create references to derived model member functions

const PtrList<volScalarField>& J = aerosol.J();
const PtrList<volScalarField>& S = aerosol.S();

const PtrList<volScalarField>& DY = diffusion.DY();
const PtrList<volScalarField>& DM = brownianDiffusion.DM();
const surfaceScalarField& phic = diffusion.phic();
const volScalarField& kEff = conductivity.kEff();

// Read velocity field

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

volScalarField DpDt
(
    IOobject
    (
        "DpDt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("DpDt", dimPressure/dimTime, 0.0)
);

volScalarField rhoStar
(
    IOobject
    (
        "rhoStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoStar", dimDensity, 0.0)
);

// Start with implicit Euler if phi doesn't exist, for consistency of phi, U and rho

IOobject phiHeader
(
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

Switch doEuler(false);

if (!phiHeader.headerOk())
{
    doEuler = true;

    Info << "Running first time step as implicit Euler." << endl;
}

// Create face flux fields

surfaceScalarField phi
(
    phiHeader,
    mesh,
    dimensionedScalar("phi", dimMass/dimTime, 0.0)
);

// Droplet fluxes

const PtrList<surfaceScalarField>& phid = aerosol.phid();

// Flux fields

surfaceScalarField phidRho
(
    IOobject
    (
        "phidRho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phidRho", phi.dimensions(), 0.0)
);

PtrList<surfaceScalarField> phiM(M.size());
PtrList<surfaceScalarField> phiY(Y.size());
PtrList<surfaceScalarField> phiZ(Z.size());

forAll(M, i)
{
    word name = M[i].name();

    phiM.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                word("phi." + name),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("phiM", phi.dimensions()/dimMass, 0.0)
        )
    );
}

forAll(species, j)
{
    word name = species.keys()[j];

    phiY.set
    (
        j,
        new surfaceScalarField
        (
            IOobject
            (
                word("phi." + name + "Y"),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("phiY", phi.dimensions(), 0.0)
        )
    );

    phiZ.set
    (
        j,
        new surfaceScalarField
        (
            IOobject
            (
                word("phi." + name + "Z"),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("phiZ", phi.dimensions(), 0.0)
        )
    );
}


// Turbulence model

autoPtr<turbModel> tTurbulence(turbModel::New(mesh, thermo));
turbModel& turbulence = tTurbulence();
const volScalarField& kTurb = turbulence.kTurb();
const volScalarField& muTurb = turbulence.muTurb();
const volScalarField& muEff = turbulence.muEff();

// Create Crank-Nicolson source terms

PtrList<volScalarField> SddtZ(thermo.nSpecies());
PtrList<volScalarField> SddtY(thermo.nSpecies());
PtrList<volScalarField> SddtYD(thermo.nSpecies());

forAll(species, j)
{
    word speciesName = species.keys()[j];

    SddtZ.set
    (
        j,
        new volScalarField
        (
            IOobject
            (
                speciesName + "ddtZ",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("ddtZ", dimDensity/dimTime, 0.0)
        )
    );

    SddtY.set
    (
        j,
        new volScalarField
        (
            IOobject
            (
                speciesName + "ddtY",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("ddtY", dimDensity/dimTime, 0.0)
        )
    );

    SddtYD.set
    (
        j,
        new volScalarField
        (
            IOobject
            (
                speciesName + "ddtYD",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("ddtYD", dimDensity/dimTime, 0.0)
        )
    );
}

PtrList<volScalarField> SddtM(M.size());

forAll(M, i)
{
    word sectionName = M[i].name();

    SddtM.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ddt" + sectionName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("ddtM", dimless/dimTime/dimVolume, 0.0)
        )
    );
}

volScalarField SddtRhoStar
(
    IOobject
    (
        "ddtRhoStar",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("ddtRhoStar", dimDensity/dimTime, 0)
);

#include "createAveragingFields.H"

// External pressure gradient

Foam::dimensionedVector externalGradP
(
    "externalGradP",
    dimPressure/dimLength,
    vector(0, 0, 0)
);

if(args.options().found("externalGradP"))
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
   );

   externalGradP = transportProperties.lookup("externalGradP");

   Info << endl << "External pressure gradient value set to: " << externalGradP << endl;
}

// Monitor field for the CFL number

volScalarField CFL
(
    IOobject
    (
        "CFL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CFL", dimless, 0.0)
);
