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

#include "semiIncompressible.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermos
{
    makeFluidThermoTypes(semiIncompressible, fluidThermo);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermos::semiIncompressible::semiIncompressible
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

Foam::fluidThermos::semiIncompressible::~semiIncompressible()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidThermos::semiIncompressible::updateRho()
{
    tmp<volScalarField> tRhoY = rhoY();
    volScalarField& rhoY = tRhoY();

    tmp<volScalarField> tRhoZ = rhoZ();
    volScalarField& rhoZ = tRhoZ();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    tmp<volScalarField> tZtot = Ztot();
    volScalarField& Ztot = tZtot();

    rhoY.max(SMALL);
    rhoZ.max(SMALL);
    Ytot.max(0.0);
    Ztot.max(0.0);

    rho_ == 1.0 / (Ytot/rhoY + Ztot/rhoZ);
}

void Foam::fluidThermos::semiIncompressible::updatePsi()
{
    psi_ *= 0;
    psi_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::semiIncompressible::rhoMix()
{
    tmp<volScalarField> tRho(new volScalarField("rhoMix", rho_));

    return tRho;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::semiIncompressible::rhoY()
{
    tmp<volScalarField> tRhoY(new volScalarField("rhoY", rho_));
    volScalarField& rhoY = tRhoY();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    Ytot.max(0.0);

    const scalarField& TCells = T_.internalField();

    PtrList<DataEntry<scalar> > dataEntriesListRho_v = getProperty("rho", VAPOR);

    volScalarField rhoInv
    (
        IOobject
        (
            "rhoInv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoInv", dimless/dimDensity, 0.0)
    );

    forAll(species_, j)
    {
        const scalarField& YCells = Y_[j].internalField();

        forAll(rhoInv.internalField(), celli)
        {
            rhoInv.internalField()[celli] +=
                max(YCells[celli], 0.0) / dataEntriesListRho_v[j].value(TCells[celli]);
        }

        forAll(rhoInv.boundaryField(), patchi)
        {
            const scalarField& TPatch = T_.boundaryField()[patchi];
            scalarField& rhoInvPatch = rhoInv.boundaryField()[patchi];
            const scalarField& YPatch = Y_[j].boundaryField()[patchi];

            forAll(rhoInvPatch, facei)
            {
                rhoInvPatch[facei] +=
                    max(YPatch[facei], 0.0) / dataEntriesListRho_v[j].value(TPatch[facei]);
            }
        }
    }

    rhoInv.max(SMALL);

    rhoY == Ytot / rhoInv;

    return tRhoY;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::semiIncompressible::rhoZ()
{
    tmp<volScalarField> tRhoZ(new volScalarField("rhoZ", rho_));
    volScalarField& rhoZ = tRhoZ();

    tmp<volScalarField> tZtot = Ztot();
    volScalarField& Ztot = tZtot();

    Ztot.max(0.0);

    const scalarField& TCells = T_.internalField();

    PtrList<DataEntry<scalar> > dataEntriesListRho_l = getProperty("rho", LIQUID);

    volScalarField rhoInv
    (
        IOobject
        (
            "rhoInv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoInv", dimless/dimDensity, 0.0)
    );

    forAll(species_, j)
    {
        const scalarField& ZCells = Z_[j].internalField();

        forAll(rhoInv.internalField(), celli)
        {
            rhoInv.internalField()[celli] +=
                max(ZCells[celli], 0.0) / dataEntriesListRho_l[j].value(TCells[celli]);
        }

        forAll(rhoInv.boundaryField(), patchi)
        {
            const scalarField& TPatch = T_.boundaryField()[patchi];
            scalarField& rhoInvPatch = rhoInv.boundaryField()[patchi];
            const scalarField& ZPatch = Z_[j].boundaryField()[patchi];

            forAll(rhoInvPatch, facei)
            {
                rhoInvPatch[facei] +=
                    max(ZPatch[facei], 0.0) / dataEntriesListRho_l[j].value(TPatch[facei]);
            }
        }
    }

    rhoInv.max(SMALL);

    rhoZ == Ztot / rhoInv;

    return tRhoZ;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::semiIncompressible::psiY()
{
    tmp<volScalarField> tPsi(new volScalarField("psiY", psi_));

    return tPsi;
}

void Foam::fluidThermos::semiIncompressible::updateCvEff()
{
    PtrList< DataEntry<scalar> >& dataEntriesListCv_v = getPropertyIfPresent("Cv", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCv_l = getPropertyIfPresent("Cv", LIQUID);

    PtrList< DataEntry<scalar> >& dataEntriesListCp_v = getPropertyIfPresent("Cp", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCp_l = getPropertyIfPresent("Cp", LIQUID);

    PtrList< DataEntry<scalar> >& dataEntriesListGamma_v = getPropertyIfPresent("gamma", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListGamma_l = getPropertyIfPresent("gamma", LIQUID);

    // Internal field

    scalarField& CvEffCells = CvEff_.internalField();
    const scalarField& TCells = T_.internalField();

    CvEffCells = 0.0;

    forAll(species_, i)
    {
        word speciesName = species_.keys()[i];

        if (propertyFound(i, "Cv", VAPOR))
        {
            forAll(CvEffCells, cellj)
            {
                CvEffCells[cellj] +=
                    max(Y_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCv_v[i].value(TCells[cellj]);
            }
        }
        else
        {
            forAll(CvEffCells, cellj)
            {
                CvEffCells[cellj] +=
                    max(Y_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCp_v[i].value(TCells[cellj]) /
                    dataEntriesListGamma_v[i].value(TCells[cellj]);
            }
        }

        if (propertyFound(i, "Cv", LIQUID))
        {
            forAll(CvEffCells, cellj)
            {
                CvEffCells[cellj] +=
                    max(Z_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCv_l[i].value(TCells[cellj]);
            }
        }
        else
        {
            forAll(CvEffCells, cellj)
            {
                CvEffCells[cellj] +=
                    max(Z_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCp_l[i].value(TCells[cellj]) /
                    dataEntriesListGamma_l[i].value(TCells[cellj]);
            }
        }
    }

    // Boundary fields

    forAll(CvEff_.boundaryField(), patchi)
    {
        scalarField& T = T_.boundaryField()[patchi];
        scalarField& CvEffPatch = CvEff_.boundaryField()[patchi];

        CvEffPatch = 0.0;

        forAll(species_, i)
        {
            word speciesName = species_.keys()[i];

            const scalarField& Y = Y_[i].boundaryField()[patchi];

            if (propertyFound(i, "Cv", VAPOR))
            {
                forAll(CvEffPatch, facej)
                {
                    CvEffPatch[facej] +=
                        max(Y[facej], 0.0) * dataEntriesListCv_v[i].value(T[facej]);
                }
            }
            else
            {
                forAll(CvEffPatch, facej)
                {
                    CvEffPatch[facej] +=
                        max(Y[facej], 0.0) * dataEntriesListCp_v[i].value(T[facej]) /
                        dataEntriesListGamma_v[i].value(T[facej]);
                }
            }

            const scalarField& Z = Z_[i].boundaryField()[patchi];

            if (propertyFound(i, "Cv", LIQUID))
            {
                forAll(CvEffPatch, facej)
                {
                    CvEffPatch[facej] +=
                        max(Z[facej], 0.0) * dataEntriesListCv_l[i].value(T[facej]);
                }
            }
            else
            {
                forAll(CvEffPatch, facej)
                {
                    CvEffPatch[facej] +=
                        max(Z[facej], 0.0) * dataEntriesListCp_l[i].value(T[facej]) /
                        dataEntriesListGamma_l[i].value(T[facej]);
                }
            }
        }
    }
}

void Foam::fluidThermos::semiIncompressible::updateCpEff()
{
    PtrList< DataEntry<scalar> >& dataEntriesListCp_v = getPropertyIfPresent("Cp", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCp_l = getPropertyIfPresent("Cp", LIQUID);

    PtrList< DataEntry<scalar> >& dataEntriesListCv_v = getPropertyIfPresent("Cv", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListCv_l = getPropertyIfPresent("Cv", LIQUID);

    PtrList< DataEntry<scalar> >& dataEntriesListGamma_v = getPropertyIfPresent("gamma", VAPOR);
    PtrList< DataEntry<scalar> >& dataEntriesListGamma_l = getPropertyIfPresent("gamma", LIQUID);

    // Internal field

    scalarField& CpEffCells = CpEff_.internalField();
    const scalarField& TCells = T_.internalField();

    CpEffCells = 0.0;

    forAll(species_, i)
    {
        word speciesName = species_.keys()[i];

        if (propertyFound(i, "Cp", VAPOR))
        {
            forAll(CpEffCells, cellj)
            {
                CpEffCells[cellj] +=
                    max(Y_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCp_v[i].value(TCells[cellj]);
            }
        }
        else
        {
            forAll(CpEffCells, cellj)
            {
                CpEffCells[cellj] +=
                    max(Y_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCv_v[i].value(TCells[cellj]) *
                    dataEntriesListGamma_v[i].value(TCells[cellj]);
            }
        }

        if (propertyFound(i, "Cp", LIQUID))
        {
            forAll(CpEffCells, cellj)
            {
                CpEffCells[cellj] +=
                    max(Z_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCp_l[i].value(TCells[cellj]);
            }
        }
        else
        {
            forAll(CpEffCells, cellj)
            {
                CpEffCells[cellj] +=
                    max(Z_[i].internalField()[cellj], 0.0) *
                    dataEntriesListCv_l[i].value(TCells[cellj]) *
                    dataEntriesListGamma_l[i].value(TCells[cellj]);
            }
        }
    }

    // Boundary fields

    forAll(CpEff_.boundaryField(), patchi)
    {
        scalarField& T = T_.boundaryField()[patchi];
        scalarField& CpEffPatch = CpEff_.boundaryField()[patchi];

        CpEffPatch = 0.0;

        forAll(species_, i)
        {
            word speciesName = species_.keys()[i];

            const scalarField& Y = Y_[i].boundaryField()[patchi];

            if (propertyFound(i, "Cp", VAPOR))
            {
                forAll(CpEffPatch, facej)
                {
                    CpEffPatch[facej] +=
                        max(Y[facej], 0.0) * dataEntriesListCp_v[i].value(T[facej]);
                }
            }
            else
            {
                forAll(CpEffPatch, facej)
                {
                    CpEffPatch[facej] +=
                        max(Y[facej], 0.0) * dataEntriesListCv_v[i].value(T[facej]) *
                        dataEntriesListGamma_v[i].value(T[facej]);
                }
            }

            const scalarField& Z = Z_[i].boundaryField()[patchi];

            if (propertyFound(i, "Cp", LIQUID))
            {
                forAll(CpEffPatch, facej)
                {
                    CpEffPatch[facej] +=
                        max(Z[facej], 0.0) * dataEntriesListCp_l[i].value(T[facej]);
                }
            }
            else
            {
                forAll(CpEffPatch, facej)
                {
                    CpEffPatch[facej] +=
                        max(Z[facej], 0.0) * dataEntriesListCv_l[i].value(T[facej]) *
                        dataEntriesListGamma_l[i].value(T[facej]);
                }
            }
        }
    }
}

Foam::tmp<Foam::scalarField>
Foam::fluidThermos::semiIncompressible::CvEff
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<Foam::scalarField> tCvEff(new scalarField(T.size()));
    scalarField& CvEff = tCvEff();

    CvEff = 1.0;

    return tCvEff;
}

Foam::tmp<Foam::scalarField> Foam::fluidThermos::semiIncompressible::CvEff
(
    const scalarField& T,
    labelList cells
) const
{
    tmp<Foam::scalarField> tCvEff(new scalarField(T.size()));
    scalarField& CvEff = tCvEff();

    CvEff = 1.0;

    return tCvEff;
}

Foam::tmp<Foam::volScalarField>
Foam::fluidThermos::semiIncompressible::CvEff
(
    const volScalarField& T
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCvEff
    (
        new volScalarField
        (
            IOobject
            (
                "CvEff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            CvEff_.dimensions()
        )
    );

    volScalarField& CvEff = tCvEff();

    // Internal

    CvEff.internalField() = this->CvEff(T.internalField());

    // Boundary patches

    forAll(CvEff.boundaryField(), patchi)
    {
        CvEff.boundaryField()[patchi] =
            this->CvEff(T.boundaryField()[patchi], patchi);
    }

    return tCvEff;
}

Foam::tmp<Foam::volScalarField>
Foam::fluidThermos::semiIncompressible::rhoLiquid()
{
    const PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        getProperty("rho", fluidThermo::LIQUID);

    tmp<volScalarField> tRho
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLiquid",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("rhoLiquid", dimDensity, 0.0)
        )
    );

    tmp<volScalarField> tZtotRho
    (
        new volScalarField
        (
            IOobject
            (
                "ZtotRho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("ZtotRho", dimless/dimDensity, 0.0)
        )
    );

    volScalarField& rho = tRho();
    volScalarField& ZtotRho = tZtotRho();

    tmp<volScalarField> tZtot = Ztot();
    volScalarField& Ztot = tZtot();

    Ztot.max(0.0);

    forAll(species_, j)
    {
        scalarField& ZtotRhoCells = ZtotRho.internalField();

        const scalarField& Zcells = Z_[j].internalField();
        const scalarField& Tcells = T_.internalField();

        forAll(ZtotRhoCells, i)
        {
            ZtotRhoCells[i] += max(Zcells[i], 0.0) / dataEntriesListRho_l[j].value(Tcells[i]);
        }

        forAll(ZtotRho.boundaryField(), patchi)
        {
            scalarField& ZtotRhoPatch = ZtotRho.boundaryField()[patchi];

            const scalarField& Zpatch = Z_[j].boundaryField()[patchi];
            const scalarField& Tpatch = T_.boundaryField()[patchi];

            forAll(ZtotRhoPatch, i)
            {
                ZtotRhoPatch[i] += max(Zpatch[i], 0.0) / dataEntriesListRho_l[j].value(Tpatch[i]);
            }
        }
    }

    rho == Ztot/stabilise(ZtotRho, dimensionedScalar("ZtotRho", dimless/dimDensity, SMALL));

    return tRho;
}

Foam::tmp<Foam::volScalarField>
Foam::fluidThermos::semiIncompressible::rhoVapor()
{
    const PtrList<DataEntry<scalar> > dataEntriesListRho_v =
        getProperty("rho", fluidThermo::VAPOR);

    tmp<volScalarField> tRho
    (
        new volScalarField
        (
            IOobject
            (
                "rhoVapor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("rhoVapor", dimDensity, 0.0)
        )
    );

    tmp<volScalarField> tYtotRho
    (
        new volScalarField
        (
            IOobject
            (
                "YtotRho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("YtotRho", dimless/dimDensity, 0.0)
        )
    );

    volScalarField& rho = tRho();
    volScalarField& YtotRho = tYtotRho();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    Ytot.max(0.0);

    forAll(species_, j)
    {
        scalarField& YtotRhoCells = YtotRho.internalField();

        const scalarField& Ycells = Y_[j].internalField();
        const scalarField& Tcells = T_.internalField();

        forAll(YtotRhoCells, i)
        {
            YtotRhoCells[i] += max(Ycells[i], 0.0) / dataEntriesListRho_v[j].value(Tcells[i]);
        }

        forAll(YtotRho.boundaryField(), patchi)
        {
            scalarField& YtotRhoPatch = YtotRho.boundaryField()[patchi];

            const scalarField& Ypatch = Y_[j].boundaryField()[patchi];
            const scalarField& Tpatch = T_.boundaryField()[patchi];

            forAll(YtotRhoPatch, i)
            {
                YtotRhoPatch[i] += max(Ypatch[i], 0.0) / dataEntriesListRho_v[j].value(Tpatch[i]);
            }
        }
    }

    rho == Ytot/stabilise(YtotRho, dimensionedScalar("YtotRho", dimless/dimDensity, SMALL));

    return tRho;
}

bool Foam::fluidThermos::semiIncompressible::read()
{
    readPropertyBoth("rho");

    return true;
}

// ************************************************************************* //
