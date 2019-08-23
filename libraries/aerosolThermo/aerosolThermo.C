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

#include "aerosolThermo.H"
#include "constants.H"
#include "harmonicPhaseMixing.H"
#include "massPhaseMixing.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aerosolThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolThermo::aerosolThermo
(
    const fvMesh& mesh
)
:
    psiThermo(mesh, word::null),
    mesh_(mesh),
    thermoCont_(nullptr),
    thermoDisp_(nullptr),
    diffusivity_(nullptr),
    phaseMix_(),
    species_(),
    activeSpecies_(),
    inactiveSpecies_(),
    activeSpeciesMap_(),
    inactiveSpeciesMap_(),
    contSpeciesMap_(),
    dispSpeciesMap_(),
    inertSpecie_(lookupType<word>("inertSpecie"))
{
    // Create specific temperature fields

    {
        volScalarField Tg
        (
            IOobject
            (
                IOobject::groupName("T", phaseNameCont()),
                mesh.time().timeName(),
                mesh
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );

        Tg.write();
    }

    {
        volScalarField Tl
        (
            IOobject
            (
                IOobject::groupName("T", phaseNameDisp()),
                mesh.time().timeName(),
                mesh
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );

        Tl.write();
    }

    // Create continuous and dispersed thermo packages

    thermoCont_ = rhoAerosolPhaseThermo::New(mesh, phaseNameCont());
    thermoDisp_ = rhoAerosolPhaseThermo::New(mesh, phaseNameDisp());

    // Set species from continuous thermo package for consistency

    wordList s(contSpecies());

    species_.transfer(s);

    wordList sActive(lookupType<wordList>("activeSpecies"));
    wordList sInactive(lookupType<wordList>("inactiveSpecies"));

    forAll(species_, j)
    {
        if (findIndex(sActive, species_[j]) != -1)
        {
            activeSpecies_.append(species_[j]);
            activeSpeciesMap_.append(j);
        }
        else
        {
            inactiveSpecies_.append(species_[j]);
            inactiveSpeciesMap_.append(j);
        }
    }

    // Set the continuous and dispersed species maps

    forAll(contSpecies(), j)
    {
        contSpeciesMap_.append(species_[contSpecies()[j]]);
    }

    forAll(dispSpecies(), j)
    {
        dispSpeciesMap_.append(species_[dispSpecies()[j]]);
    }

    // Set active and inactive species explicitly

    forAll(contSpecies(), j)
    {
        const word speciesName(contSpecies()[j]);

        if (activeSpecies_.found(speciesName))
        {
            thermoCont_->composition().setActive(j);
        }
        else
        {
            thermoCont_->composition().setInactive(j);
        }
    }

    forAll(dispSpecies(), j)
    {
        const word speciesName(dispSpecies()[j]);

        if (activeSpecies_.found(speciesName))
        {
            thermoDisp_->composition().setActive(j);
        }
        else
        {
            thermoDisp_->composition().setInactive(j);
        }
    }

    phaseMix_.set(new phaseMixing(*this));

    correct();

    diffusivity_.set(new mixtureDiffusivityModel(*this));

    // Print species arays

    Info<< "species =" << endl << "(" << endl;

    forAll(this->species(), j)
    {
        Info<< token::TAB << this->species()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "contSpecies =" << endl << "(" << endl;

    forAll(this->contSpecies(), j)
    {
        Info<< token::TAB << this->contSpecies()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "dispSpecies =" << endl << "(" << endl;

    forAll(this->dispSpecies(), j)
    {
        Info<< token::TAB << this->dispSpecies()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "activeSpecies =" << endl << "(" << endl;

    forAll(this->activeSpecies(), j)
    {
        Info<< token::TAB << this->activeSpecies()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "inactiveSpecies =" << endl << "(" << endl;

    forAll(this->inactiveSpecies(), j)
    {
        Info<< token::TAB << this->inactiveSpecies()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "activeSpeciesMap =" << endl << "(" << endl;

    forAll(this->activeSpeciesMap(), j)
    {
        Info<< token::TAB << this->activeSpeciesMap()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "inactiveSpeciesMap =" << endl << "(" << endl;

    forAll(this->inactiveSpeciesMap(), j)
    {
        Info<< token::TAB << this->inactiveSpeciesMap()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "contSpeciesMap =" << endl << "(" << endl;

    forAll(this->contSpeciesMap(), j)
    {
        Info<< token::TAB << this->contSpeciesMap()[j] << endl;
    }

    Info<< ")" << endl;

    Info<< "dispSpeciesMap =" << endl << "(" << endl;

    forAll(this->dispSpeciesMap(), j)
    {
        Info<< token::TAB << this->dispSpeciesMap()[j] << endl;
    }

    Info<< ")" << endl << endl;

    // Set multivariate fields

    forAll(this->contSpecies(), j)
    {
        fields_.add(this->Y()[j]);
        fieldsY_.add(this->Y()[j]);
    }

    forAll(this->dispSpecies(), j)
    {
        fields_.add(this->Z()[j]);
        fieldsZ_.add(this->Z()[j]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolThermo::~aerosolThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolThermo::correctThermo()
{
    thermoCont_->T() = T_;
    thermoCont_->he() = thermoCont_->he(p_, T_);
    thermoCont_->correct();

    thermoDisp_->T() = T_;
    thermoDisp_->he() = thermoDisp_->he(p_, T_);
    thermoDisp_->correct();
}


void Foam::aerosolThermo::correct()
{
    const volScalarField rho
    (
        harmonicPhaseMixing(*this).mix
        (
            thermoCont_->rho(),
            thermoDisp_->rho()
        )
    );

    psi_ = rho/p_;

    mu_ = phaseMix_->viscosity().mix
    (
        thermoCont_->mu(),
        thermoDisp_->mu()
    );

    alpha_ = phaseMix_->conductivity().mix
    (
        thermoCont_->alpha(),
        thermoDisp_->alpha()
    );
}


bool Foam::aerosolThermo::incompressible() const
{
    return thermoCont_->incompressible() && thermoDisp_->incompressible();
}


bool Foam::aerosolThermo::isochoric() const
{
    return thermoCont_->isochoric() && thermoDisp_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return massPhaseMixing(*this).mix
    (
        thermoCont_->he(p, T),
        thermoDisp_->he(p, T)
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return massPhaseMixing(*this).mix
    (
        thermoCont_->he(p, T, cells),
        thermoDisp_->he(p, T, cells),
        cells
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return massPhaseMixing(*this).mix
    (
        thermoCont_->he(p, T, patchi),
        thermoDisp_->he(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::hc() const
{
    return massPhaseMixing(*this).mix
    (
        thermoCont_->hc(),
        thermoDisp_->hc()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::rho() const
{
    return psiThermo::rho();
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::rho
(
    const label patchi
) const
{
    return psiThermo::rho(patchi);
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::rho
(
    const labelList& cells
) const
{
    tmp<scalarField> trho(new scalarField(cells.size(), 0.0));

    scalarField& rho = trho.ref();

    forAll(cells, i)
    {
        rho[i] = p_[cells[i]]*psi_[cells[i]];
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::Cp() const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cp(),
        thermoDisp_->Cp()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cp(p, T, patchi),
        thermoDisp_->Cp(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::Cv() const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cv(),
        thermoDisp_->Cv()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cv(p, T, patchi),
        thermoDisp_->Cv(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::gamma() const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->gamma(),
        thermoDisp_->gamma()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->gamma(p, T, patchi),
        thermoDisp_->gamma(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::lambda() const
{
    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar kB(constant::physicoChemical::k);
    const dimensionedScalar NA(constant::physicoChemical::NA);

    const volScalarField m(thermoCont_->WMix()*0.001/NA);

    return thermoCont_->mu()/p_*Foam::sqrt(pi*kB*T_/(2*m));
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::Cpv() const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cpv(),
        thermoDisp_->Cpv()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->Cpv(p, T, patchi),
        thermoDisp_->Cpv(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::CpByCpv() const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->CpByCpv(),
        thermoDisp_->CpByCpv()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return phaseMix_->heatCapacity().mix
    (
        thermoCont_->CpByCpv(p, T, patchi),
        thermoDisp_->CpByCpv(p, T, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::nu
(
    const label patchi
) const
{
    return mu(patchi)/rho(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::kappa() const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->kappa(),
        thermoDisp_->kappa()
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::kappa
(
    const label patchi
) const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->kappa(patchi),
        thermoDisp_->kappa(patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->kappaEff(alphat),
        thermoDisp_->kappaEff(alphat)
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->kappaEff(alphat, patchi),
        thermoDisp_->kappaEff(alphat, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->alphaEff(alphat),
        thermoDisp_->alphaEff(alphat)
    );
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return phaseMix_->conductivity().mix
    (
        thermoCont_->alphaEff(alphat, patchi),
        thermoDisp_->alphaEff(alphat, patchi),
        patchi
    );
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::sumY() const
{
    tmp<volScalarField> tsumY
    (
        new volScalarField
        (
            IOobject
            (
                "sumY",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Y()[0]*0.0 // Necessary to have proper calculated patches!
        )
    );

    volScalarField& sumY = tsumY.ref();

    forAll(Y(), j)
    {
        sumY += Y()[j];
    }

    return tsumY;
}

Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::sumZ() const
{
    tmp<volScalarField> tsumZ
    (
        new volScalarField
        (
            IOobject
            (
                "sumZ",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Z()[0]*0.0 // Necessary to have proper calculated patches!
        )
    );

    volScalarField& sumZ = tsumZ.ref();

    forAll(Z(), j)
    {
        sumZ += Z()[j];
    }

    return tsumZ;
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::sumY
(
    const label patchi
) const
{
    tmp<scalarField> tsumY(new scalarField(Y()[0].boundaryField()[patchi]));

    if (Y().size() > 1)
    {
        scalarField& sumY = tsumY.ref();

        for (label j = 1; j < Y().size(); j++)
        {
            sumY += Y()[j].boundaryField()[patchi];
        }
    }

    return tsumY;
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::sumZ
(
    const label patchi
) const
{
    tmp<scalarField> tsumZ(new scalarField(Z()[0].boundaryField()[patchi]));

    if (Z().size() > 1)
    {
        scalarField& sumZ = tsumZ.ref();

        for (label j = 1; j < Z().size(); j++)
        {
            sumZ += Z()[j].boundaryField()[patchi];
        }
    }

    return tsumZ;
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::sumY
(
    const labelList& cells
) const
{
    tmp<scalarField> tsumY(new scalarField(cells.size(), 0.0));

    scalarField& sumY = tsumY.ref();

    forAll(Y(), j)
    {
        forAll(cells, i)
        {
            sumY[i] += Y()[j].field()[cells[i]];
        }
    }

    return tsumY;
}


Foam::tmp<Foam::scalarField> Foam::aerosolThermo::sumZ
(
    const labelList& cells
) const
{
    tmp<scalarField> tsumZ(new scalarField(cells.size(), 0.0));

    scalarField& sumZ = tsumZ.ref();

    forAll(Z(), j)
    {
        forAll(cells, i)
        {
            sumZ[i] += Z()[j].field()[cells[i]];
        }
    }

    return tsumZ;
}


Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::W() const
{
    return this->WMix();
}

Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::alphahe() const
{
    NotImplemented;
    return tmp<volScalarField>();
}

Foam::tmp<Foam::scalarField> Foam::aerosolThermo::alphahe
(
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>();
}

Foam::tmp<Foam::volScalarField> Foam::aerosolThermo::WMix() const
{
    return harmonicPhaseMixing(*this).mix
    (
        thermoCont_->WMix(),
        thermoDisp_->WMix()
    );
}

Foam::tmp<Foam::scalarField> Foam::aerosolThermo::WMix
(
    const label patchi
) const
{
    return harmonicPhaseMixing(*this).mix
    (
        thermoCont_->WMix(patchi),
        thermoDisp_->WMix(patchi),
        patchi
    );
}

Foam::tmp<Foam::scalarField> Foam::aerosolThermo::WMix
(
    const labelList& cells
) const
{
    return harmonicPhaseMixing(*this).mix
    (
        thermoCont_->WMix(cells),
        thermoDisp_->WMix(cells),
        cells
    );
}


Foam::PtrList<Foam::scalarField> Foam::aerosolThermo::pSat
(
    const speciesTable& species
)
{
    PtrList<scalarField> pSat(species.size());

    const scalarField& T = T_;

    forAll(species, j)
    {
        pSat.set
        (
            j,
            new scalarField
            (
                thermoCont_->property(species[j], "pSat").value(T)
            )
        );
    }

    return pSat;
}


Foam::PtrList<Foam::scalarField> Foam::aerosolThermo::rhoDisp
(
    const speciesTable& species
) const
{
    PtrList<scalarField> rhoDisp(species.size());

    const scalarField& T = T_;
    const scalarField& p = p_;

    const basicSpecieMixture& comp = thermoDisp_->composition();

    forAll(species, j)
    {
        rhoDisp.set(j, new scalarField(T.size(), 0.0));

        forAll(T_, celli)
        {
            rhoDisp[j][celli] = comp.rho(j, p[celli], T[celli]);
        }
    }

    return rhoDisp;
}


Foam::PtrList<Foam::scalarField> Foam::aerosolThermo::rhoCont
(
    const speciesTable& species
) const
{
    PtrList<scalarField> rhoCont(species.size());

    const scalarField& T = T_;
    const scalarField& p = p_;

    const basicSpecieMixture& comp = thermoCont_->composition();

    forAll(species, j)
    {
        rhoCont.set(j, new scalarField(T.size(), 0.0));

        forAll(T_, celli)
        {
            rhoCont[j][celli] = comp.rho(j, p[celli], T[celli]);
        }
    }

    return rhoCont;
}


Foam::PtrList<Foam::scalarField> Foam::aerosolThermo::sigma
(
    const speciesTable& species
)
{
    PtrList<scalarField> sigma(species.size());

    const scalarField& T = T_;

    forAll(species, j)
    {
        sigma.set
        (
            j,
            new scalarField
            (
                thermoDisp_->property(species[j], "sigma").value(T)
            )
        );
    }

    return sigma;
}

bool Foam::aerosolThermo::read()
{
    if (psiThermo::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
