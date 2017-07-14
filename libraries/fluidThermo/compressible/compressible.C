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

#include "compressible.H"
#include "makeFluidThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermos
{
    makeFluidThermoTypes(compressible, fluidThermo);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermos::compressible::compressible
(
    const fvMesh& mesh
)
:
    fluidThermo(mesh),
    mesh_(mesh)
{
    this->read();

    // Compute fields

    if (!psi_.headerOk()) updatePsi();
    if (!rho_.headerOk()) updateRho();

    // Necessary for resuming a simulation, when updatePsi() is called before pEqn:
    psi_.oldTime();

    updateCvEff();
    updateCpEff();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermos::compressible::~compressible()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidThermos::compressible::updateRho()
{
    rho_ = psi_ * (p1_ + p0_);
}

void Foam::fluidThermos::compressible::updatePsi()
{
    tmp<volScalarField> tPsiY = psiY();
    volScalarField& psiY = tPsiY();

    tmp<volScalarField> tRhoZ = rhoZ();
    volScalarField& rhoZ = tRhoZ();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    tmp<volScalarField> tZtot = Ztot();
    volScalarField& Ztot = tZtot();

    psiY.max(SMALL);
    rhoZ.max(SMALL);
    Ytot.max(0.0);
    Ztot.max(0.0);

    psi_ == 1.0 / ((Ytot / psiY) + ((p1_ + p0_) * Ztot / rhoZ));
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::rhoMix()
{
    tmp<volScalarField> tRho(new volScalarField("rhoMix", rho_));
    volScalarField& rho = tRho();

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

    rho == 1.0 / (Ytot/rhoY + Ztot/rhoZ);

    return tRho;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::rhoY()
{
    tmp<volScalarField> tRhoY(new volScalarField("rhoY", rho_));
    volScalarField& rhoY = tRhoY();

    rhoY == psiY() * (p1_ + p0_);

    return tRhoY;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::rhoZ()
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

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::psiY()
{
    tmp<volScalarField> tPsi
    (
        new volScalarField
        (
            IOobject
            (
                "psiY",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("psiY", psi_.dimensions(), 0.0)
        )
    );

    volScalarField& psi = tPsi();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    Ytot.max(0.0);

    const scalarField& TCells = T_.internalField();

    PtrList< DataEntry<scalar> >& dataEntriesListPsi_v = getProperty("psi", VAPOR);

    volScalarField psiInv
    (
        IOobject
        (
            "psiInv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("psiInv", dimPressure/dimDensity, 0.0)
    );

    forAll(species_, j)
    {
        const scalarField& YCells = Y_[j].internalField();

        forAll(psiInv.internalField(), celli)
        {
            psiInv.internalField()[celli] +=
                max(YCells[celli], 0.0) / dataEntriesListPsi_v[j].value(TCells[celli]);
        }

        forAll(psiInv.boundaryField(), patchi)
        {
            const scalarField& TPatch = T_.boundaryField()[patchi];
            scalarField& psiInvPatch = psiInv.boundaryField()[patchi];
            const scalarField& YPatch = Y_[j].boundaryField()[patchi];

            forAll(psiInvPatch, facei)
            {
                psiInvPatch[facei] +=
                    max(YPatch[facei], 0.0) / dataEntriesListPsi_v[j].value(TPatch[facei]);
            }
        }
    }

    psiInv.max(SMALL);

    psi == Ytot / psiInv;

    return tPsi;
}

void Foam::fluidThermos::compressible::updateCvEff()
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

void Foam::fluidThermos::compressible::updateCpEff()
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

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::rhoLiquid()
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

Foam::tmp<Foam::volScalarField> Foam::fluidThermos::compressible::rhoVapor()
{
    const PtrList<DataEntry<scalar> > dataEntriesListPsi_v =
        getProperty("psi", fluidThermo::VAPOR);

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
        const scalarField& p1cells = p1_.internalField();

        forAll(YtotRhoCells, i)
        {
            YtotRhoCells[i] += max(Ycells[i], 0.0)
                             / (p1cells[i] + p0_.value())
                             / dataEntriesListPsi_v[j].value(Tcells[i]);
        }

        forAll(YtotRho.boundaryField(), patchi)
        {
            scalarField& YtotRhoPatch = YtotRho.boundaryField()[patchi];

            const scalarField& Ypatch = Y_[j].boundaryField()[patchi];
            const scalarField& Tpatch = T_.boundaryField()[patchi];
            const scalarField& p1patch = p1_.boundaryField()[patchi];

            forAll(YtotRhoPatch, i)
            {
                YtotRhoPatch[i] += max(Ypatch[i], 0.0)
                                 / (p1patch[i] + p0_.value())
                                 / dataEntriesListPsi_v[j].value(Tpatch[i]);
            }
        }
    }

    rho == Ytot/stabilise(YtotRho, dimensionedScalar("YtotRho", dimless/dimDensity, SMALL));

    return tRho;
}

bool Foam::fluidThermos::compressible::read()
{
    // Vapors are assumed to be compressible according to some psi. The density
    // is computed with rho = psi*p.

    readProperty("psi", VAPOR);

    // The density of liquid is specified with rho. Psi is computed with
    // psi = rho/p.

    readProperty("rho", LIQUID);

    return true;
}

// ************************************************************************* //
