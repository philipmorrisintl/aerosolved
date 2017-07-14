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

#include "incompressible.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermos
{
    makeFluidThermoTypes(incompressible, fluidThermo);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermos::incompressible::incompressible
(
    const fvMesh& mesh
)
:
    fluidThermo(mesh),
    mesh_(mesh)
{
    read();

    // Compute fields

    updatePsi();
    updateRho();
    updateCvEff();
    updateCpEff();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermos::incompressible::~incompressible()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidThermos::incompressible::updateRho()
{
    rho_ == dimensionedScalar("rho", dimDensity, rhoValue_);
    rho_.correctBoundaryConditions();
}

void Foam::fluidThermos::incompressible::updatePsi()
{
    psi_ *= 0;
    psi_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::incompressible::rhoMix()
{
    tmp<volScalarField> tRho(new volScalarField("rhoMix", rho_));

    return tRho;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::incompressible::rhoY()
{
    tmp<volScalarField> tRhoY(new volScalarField("rhoY", rho_));

    return tRhoY;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::incompressible::rhoZ()
{
    tmp<volScalarField> tRhoZ(new volScalarField("rhoZ", rho_));

    volScalarField& rhoZ = tRhoZ();

    rhoZ *= 0.0;
    rhoZ.correctBoundaryConditions();

    return tRhoZ;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::incompressible::psiY()
{
    tmp<volScalarField> tPsi(new volScalarField("psiY", psi_));

    return tPsi;
}

void Foam::fluidThermos::incompressible::updateCvEff()
{
    PtrList< DataEntry<scalar> >& dataEntriesListCv_v = getPropertyIfPresent("Cv", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCp_v = getPropertyIfPresent("Cp", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListGamma_v = getPropertyIfPresent("gamma", VAPOR);

    // Internal field

    scalarField& CvEffCells = CvEff_.internalField();
    const scalarField& TCells = T_.internalField();

    if (propertyFound(0, "Cv", VAPOR))
    {
        forAll(CvEffCells, cellj)
        {
            CvEffCells[cellj] = dataEntriesListCv_v[0].value(TCells[cellj]);
        }
    }
    else
    {
        forAll(CvEffCells, cellj)
        {
            CvEffCells[cellj] =
                dataEntriesListCp_v[0].value(TCells[cellj]) /
                dataEntriesListGamma_v[0].value(TCells[cellj]);
        }
    }

    // Boundary fields

    forAll(CvEff_.boundaryField(), patchi)
    {
        scalarField& T = T_.boundaryField()[patchi];
        scalarField& CvEffPatch = CvEff_.boundaryField()[patchi];

        if (propertyFound(0, "Cv", VAPOR))
        {
            forAll(CvEffPatch, facej)
            {
                CvEffPatch[facej] = dataEntriesListCv_v[0].value(T[facej]);
            }
        }
        else
        {
            forAll(CvEffPatch, facej)
            {
                CvEffPatch[facej] =
                        dataEntriesListCp_v[0].value(T[facej]) /
                        dataEntriesListGamma_v[0].value(T[facej]);
            }
        }
    }
}

void Foam::fluidThermos::incompressible::updateCpEff()
{
    PtrList< DataEntry<scalar> >& dataEntriesListCp_v = getPropertyIfPresent("Cp", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCv_v = getPropertyIfPresent("Cv", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListGamma_v = getPropertyIfPresent("gamma", VAPOR);

    // Internal field

    scalarField& CpEffCells = CpEff_.internalField();
    const scalarField& TCells = T_.internalField();

    if (propertyFound(0, "Cp", VAPOR))
    {
        forAll(CpEffCells, cellj)
        {
            CpEffCells[cellj] = dataEntriesListCp_v[0].value(TCells[cellj]);
        }
    }
    else
    {
        forAll(CpEffCells, cellj)
        {
            CpEffCells[cellj] =
                    dataEntriesListCv_v[0].value(TCells[cellj]) *
                    dataEntriesListGamma_v[0].value(TCells[cellj]);
        }
    }

    // Boundary fields

    forAll(CpEff_.boundaryField(), patchi)
    {
        scalarField& T = T_.boundaryField()[patchi];
        scalarField& CpEffPatch = CpEff_.boundaryField()[patchi];

        if (propertyFound(0, "Cp", VAPOR))
        {
            forAll(CpEffPatch, facej)
            {
                CpEffPatch[facej] += dataEntriesListCp_v[0].value(T[facej]);
            }
        }
        else
        {
            forAll(CpEffPatch, facej)
            {
                CpEffPatch[facej] =
                    dataEntriesListCv_v[0].value(T[facej]) *
                    dataEntriesListGamma_v[0].value(T[facej]);
            }
        }
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fluidThermos::incompressible::rhoLiquid()
{
    tmp<volScalarField> tRho(new volScalarField("rho", rho_));

    volScalarField& rho = tRho();

    rho *= 0.0;
    rho.correctBoundaryConditions();

    return tRho;
}

Foam::tmp<Foam::volScalarField>
Foam::fluidThermos::incompressible::rhoVapor()
{
    tmp<volScalarField> tRho(new volScalarField("rho", rho_));

    return tRho;
}

bool Foam::fluidThermos::incompressible::read()
{
    params_.lookup("rhoValue") >> rhoValue_;

    return true;
}

// ************************************************************************* //
