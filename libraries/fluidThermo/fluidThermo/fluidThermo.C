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

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, dictionary)
}

const Foam::word Foam::fluidThermo::dictName("fluidThermoProperties");


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

bool Foam::fluidThermo::readPropertyThisSpecies
(
    PtrList< DataEntry<scalar> >& dataEntrieList,
    const label i,
    const word propertyName,
    const phaseType phase,
    const dictionary& speciesDict
)
{
    word speciesName = speciesDict.keys()[i];
    dictionary sDict(speciesDict.subDict(speciesName));
    word dictName = (phase == VAPOR ? "vaporProperties" : "liquidProperties");

    if (sDict.subDict(dictName).found(propertyName))
    {
        dataEntrieList.set(i, DataEntry<scalar>::New(propertyName, sDict.subDict(dictName)));
        return true;
    }

    return false;
}

void Foam::fluidThermo::readProperty
(
    const word propertyName,
    const phaseType phase,
    const dictionary& speciesDict,
    const bool force
)
{
    HashTable< PtrList<DataEntry<scalar> > > *dataEntriesPtr(NULL);
    HashTable< List<Switch> > *propertiesAvailablePtr(NULL);

    if (phase == VAPOR)
    {
        dataEntriesPtr = &dataEntriesVapor_;
        propertiesAvailablePtr = &propertiesAvailableVapor_;
    }
    else if (phase == LIQUID)
    {
        dataEntriesPtr = &dataEntriesLiquid_;
        propertiesAvailablePtr = &propertiesAvailableLiquid_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::fluidThermo::readProperty(const word propertyName, const word phase, const dictionary& speciesDict, const bool force)"
        )   << "No valid phase specified. Must be either vapor or liquid." << nl
            << exit(FatalError);
    }


    // If this entry is not yet created, do so

    if (!dataEntriesPtr->found(propertyName))
    {
        dictionary dummyDict;
        dummyDict.add("zero", "constant 0.0");

        PtrList< DataEntry<scalar> > emptyDataEntriesList(nSpecies());

        forAll(emptyDataEntriesList, i)
        {
            emptyDataEntriesList.set(i, DataEntry<scalar>::New("zero", dummyDict));
        }

        dataEntriesPtr->insert(propertyName, emptyDataEntriesList);
    }

    if (!propertiesAvailablePtr->found(propertyName))
    {
        List<Switch> emptyPropertiesAvailableList(species_.size(), false);
        propertiesAvailablePtr->insert(propertyName, emptyPropertiesAvailableList);
    }

    // Create reference into dataEntriesPtr and propertiesAvailablePtr

    PtrList<DataEntry<scalar> >& dataEntriesList = dataEntriesPtr->find(propertyName)();
    List<Switch>& propertiesAvailableList = propertiesAvailablePtr->find(propertyName)();

    // Read property

    forAll(speciesDict, i)
    {
        word speciesName = speciesDict.keys()[i];

        propertiesAvailableList[i] =
            readPropertyThisSpecies(dataEntriesList, i, propertyName, phase, speciesDict);

        if (force && !propertiesAvailableList[i])
        {
            FatalErrorIn
            (
                "Foam::fluidThermo::readProperty(const word propertyName, const word phase, const dictionary& speciesDict, const bool force)"
            )   << "Property " << propertyName << " for species " << speciesName
                << " (" << (phase == VAPOR ? "vapor" : "liquid")
                << ") is required, but was not found." << nl
                << exit(FatalError);
        }
    }
}

Foam::PtrList<Foam::DataEntry<Foam::scalar> >& Foam::fluidThermo::getProperty
(
    const word propertyName,
    const phaseType phase,
    const bool force
)
{
    HashTable< PtrList<DataEntry<scalar> > > *dataEntriesPtr(NULL);

    if (phase == VAPOR)
    {
        dataEntriesPtr = &dataEntriesVapor_;
    }
    else if (phase == LIQUID)
    {
        dataEntriesPtr = &dataEntriesLiquid_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::fluidThermo::getProperty(const word propertyName, const phaseType phase, const bool force)"
        )   << "No valid phase specified. Must be either vapor or liquid." << nl
            << exit(FatalError);
    }

    if(!dataEntriesPtr->found(propertyName))
    {
        if (force)
        {
            FatalErrorIn
            (
                "Foam::fluidThermo::getProperty(const word propertyName, const phaseType phase, const bool force)"
            )   << "Cannot get property " << propertyName
                << " (" << (phase == VAPOR ? "vapor" : "liquid") << "). You must read it first "
                << "with readProperty()." << exit(FatalError);
        }
        else
        {
            dictionary dummyDict;
            dummyDict.add("zero", "constant 0.0");

            PtrList< DataEntry<scalar> > emptyDataEntriesList(nSpecies());

            forAll(emptyDataEntriesList, i)
            {
                emptyDataEntriesList.set(i, DataEntry<scalar>::New("zero", dummyDict));
            }

            dataEntriesPtr->insert(propertyName, emptyDataEntriesList);
        }
    }

    return dataEntriesPtr->find(propertyName)();
}

void Foam::fluidThermo::getDiffusivityModels()
{
    diffusivityModels_.clear();

    diffusivityModels_.setSize(nSpecies()*(nSpecies()-1));

    label k = 0;

    forAll(species_, i)
    {
        forAll(species_, j)
        {
            if (i < j)
            {
                word a = species_.keys()[i];
                word b = species_.keys()[j];

                dictionary av = species_.subDict(a).subDict("vaporProperties");
                dictionary bv = species_.subDict(b).subDict("vaporProperties");

                if
                (
                    av.found("diffusivity") &&
                    av.subDict("diffusivity").found(b)
                )
                {
                    diffusivityModels_.set
                    (
                        k,
                        diffusivityModel::New
                        (
                            b,
                            av.subDict("diffusivity"),
                            species_,
                            i,
                            j
                        )
                    );
                }
                else if
                (
                    av.found("diffusivity") &&
                    av.subDict("diffusivity").found("D12")
                )
                {
                    diffusivityModels_.set
                    (
                        k,
                        diffusivityModel::New
                        (
                            "D12",
                            av.subDict("diffusivity"),
                            species_,
                            i,
                            j
                        )
                    );
                }
                else if
                (
                    bv.found("diffusivity") &&
                    bv.subDict("diffusivity").found(a)
                )
                {
                    diffusivityModels_.set
                    (
                        k,
                        diffusivityModel::New
                        (
                            a,
                            bv.subDict("diffusivity"),
                            species_,
                            i,
                            j
                        )
                    );
                }
                else
                {
                    FatalErrorIn
                    (
                        "Foam::fluidThermo::getDiffusivityModels()"
                    )   << "Cannot find diffusivity for " << a << " and " << b << nl
                        << exit(FatalError);
                }

                diffIndex_[i][j] = k;
                diffIndex_[j][i] = k;

                k++;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),

    // Fields to read
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p1_
    (
        IOobject
        (
            "p1",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p0_("p0", dimless, 0.0),


    // Fields to compute
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity, 0.0)
    ),
    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity/dimPressure, 0.0)
    ),
    CvEff_
    (
        IOobject
        (
            "CvEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTemperature/dimMass, 0.0)
    ),
    CpEff_
    (
        IOobject
        (
            "CpEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTemperature/dimMass, 0.0)
    ),
    W_(subDict("Species").size()),
    X_(subDict("Species").size()),
    Y_(subDict("Species").size()),
    Z_(subDict("Species").size()),
    R_(dimensionedScalar("R", dimensionSet(1,2,-2,-1,-1,0,0), 8.3144621)),
    M_(subDict("Species").size()),
    Tc_(subDict("Species").size()),
    phaseChange_(subDict("Species").size()),
    species_(subDict("Species")),
    speciesPhaseChange_(subDict("Species")),
    params_(subDict("fluidThermoModelParameters")),
    dataEntriesVapor_(0),
    dataEntriesLiquid_(0),
    propertiesAvailableVapor_(0),
    propertiesAvailableLiquid_(0),
    diffusivityModels_(0),
    diffIndex_(subDict("Species").size(), subDict("Species").size(), -1)
{
    word speciesName;

    // Order by phase change switch

    int j = 0;

    for(int i = 0; i < species_.size(); i++)
    {
        speciesName = species_.keys()[j];

        if(!species_.subDict(speciesName).lookupOrDefault<Switch>("phaseChange", true))
        {
            dictionary thisSpecies = species_.subOrEmptyDict(speciesName);

            if(species_.remove(speciesName))
            {
                species_.set(speciesName, thisSpecies);
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::fluidThermo::fluidThermo(const fvMesh& mesh)"
                )   << "cannot move species in species dictionary" << nl
                    << exit(FatalError);
            }
        }
        else
        {
            j++;
        }
    }

    // Set phase change switch and fill phase change dictionary

    speciesPhaseChange_ = species_;

    forAll(species_, i)
    {
        speciesName = species_.keys()[i];
        phaseChange_[i] = species_.subDict(speciesName).lookupOrDefault<Switch>("phaseChange", true);

        if (!phaseChange_[i])
        {
            if(!speciesPhaseChange_.remove(speciesName))
            {
                FatalErrorIn
                (
                    "Foam::fluidThermo::fluidThermo(const fvMesh& mesh)"
                )   << "cannot remove species from species phase change dictionary" << nl
                    << exit(FatalError);
            }
        }
    }

    // Read info for each species

    forAll(species_, i)
    {
        speciesName = species_.keys()[i];

        // Get vapor mass fraction
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesName + "Y",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // Create vapor mole fraction field
        W_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesName + "W",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                scalar(0.0)
            )
        );

        // Get liquid mass fraction
        Z_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesName + "Z",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // Create liquid mole fraction field
        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesName + "X",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                scalar(0.0)
            )
        );

        // Get molar mass and diffusion volume
        species_.subDict(speciesName).lookup("moleWeight") >> M_[i];

        // Get critical temperature
        species_.subDict(speciesName).lookup("Tc") >> Tc_[i];
    }

    updateMoleFractions();

    this->read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fluidThermo::read()
{
    massConservationTolerance_ =
        params_.lookupOrDefault<scalar>("massConservationTolerance", 0.1);
    massConservationRelaxation_ =
        params_.lookupOrDefault<scalar>("massConservationRelaxation", 0.1);

    params_.lookup("p0") >> p0_;

    params_.lookup("rholMin") >> rholMin_;
    params_.lookup("rholMax") >> rholMax_;
    params_.lookup("rhovMin") >> rhovMin_;
    params_.lookup("rhovMax") >> rhovMax_;

    // Read Cp, Cv and gamma if present

    readPropertyBothIfPresent("Cp");
    readPropertyBothIfPresent("Cv");
    readPropertyBothIfPresent("gamma");

    // Check if sufficient information is available. Fill the gammaPow field
    // so that Cp can be computed from Cv, and vice versa

    for(int j = 0; j <= 1; j++)
    {
        phaseType phase = phaseType(j);

        forAll(species_, i)
        {
            word speciesName = species_.keys()[i];

            if (!((propertyFound(i, "Cp", phase) && propertyFound(i, "Cv", phase)) ||
                  (propertyFound(i, "Cp", phase) && propertyFound(i, "gamma", phase)) ||
                  (propertyFound(i, "Cv", phase) && propertyFound(i, "gamma", phase))))
            {
                FatalErrorIn
                (
                    "Foam::fluidThermos::compressible::read()"
                )   << "No Cp and Cv, Cp and gamma or Cv and gamma found for " << speciesName
                    << " (" << (phase == VAPOR ? "vapor" : "liquid") << ")."
                    << exit(FatalError);
            }
        }
    }

    return true;
}

const Foam::fluidThermo& Foam::fluidThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    const fluidThermo& thermo = pf.db().lookupObject<fluidThermo>(dictName);
    return thermo;
}

void Foam::fluidThermo::rescale
(
    PtrList<volScalarField>& Y,
    PtrList<volScalarField>& Z,
    PtrList<volScalarField>& M,
    bool printInfo,
    bool clipToZero
)
{
    tmp<volScalarField> tYt
    (
        new volScalarField
        (
            IOobject
            (
                "Yt",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            scalar(0.0)
        )
    );

    volScalarField& Yt = tYt();

    // Clip the solutions to zero

    if (clipToZero)
    {
        forAll(Y, j)
        {
            Y[j].max(dimensionedScalar("zero", dimless, 0.0));
            Z[j].max(dimensionedScalar("zero", dimless, 0.0));

            Y[j].correctBoundaryConditions();
            Z[j].correctBoundaryConditions();
        }

        forAll(M, i)
        {
            M[i].max(dimensionedScalar("zero", M[i].dimensions(), 0.0));

            M[i].correctBoundaryConditions();
        }
    }

    // Compute total mass fraction

    forAll(Y, j)
    {
        Yt += (Y[j] + Z[j]);
    }

    if (debug)
    {
        Pout << "Max Yt: " << max(Yt.internalField()) << endl;
        Pout << "Min Yt: " << min(Yt.internalField()) << endl;
    }

    const scalar maxInternal(gMax(Yt));
    const scalar minInternal(gMin(Yt));

    scalar maxBoundary(1.0);
    scalar minBoundary(1.0);

    forAll(Yt.boundaryField(), patchi)
    {
        scalarField& YtPatch = Yt.boundaryField()[patchi];

        maxBoundary = max(maxBoundary, max(YtPatch));
        minBoundary = min(minBoundary, min(YtPatch));
    }

    reduce(maxBoundary, maxOp<scalar>());
    reduce(minBoundary, minOp<scalar>());

    if (printInfo)
    {
        Info << "Scaling species with relaxation = "
             << massConservationRelaxation_ << ", max(Yt) = " << max(maxInternal, maxBoundary)
             << ", min(Yt) = " << min(minInternal, minBoundary) << endl;
    }

    if (max(maxInternal, maxBoundary) > (1.0 + massConservationTolerance_))
    {
        FatalErrorIn
        (
            "Foam::fluidThermo::rescale(...)"
        )   << "Sum of the mass fractions became larger than 1+"
            << massConservationTolerance_ << ". Exiting."
            << exit(FatalError);
    }

    if (min(minInternal, minBoundary) < (1.0 - massConservationTolerance_))
    {
        FatalErrorIn
        (
            "Foam::fluidThermo::rescale(...)"
        )   << "Sum of the mass fractions became smaller than 1-"
            << massConservationTolerance_ << ". Exiting."
            << exit(FatalError);
    }

    if (massConservationRelaxation_ > SMALL)
    {
        dimensionedScalar one("one", dimless, 1.0);
        dimensionedScalar r("r", dimless, massConservationRelaxation_);

        forAll(Y, j)
        {
            Y[j] /= (one+(Yt-one)*r);
            Z[j] /= (one+(Yt-one)*r);

            Y[j].correctBoundaryConditions();
            Z[j].correctBoundaryConditions();
        }

        forAll(M, i)
        {
            M[i] /= (one+(Yt-one)*r);

            M[i].correctBoundaryConditions();
        }
    }
}

void Foam::fluidThermo::updateMoleFractions()
{
    tmp<volScalarField> tMinv
    (
        new volScalarField
        (
            IOobject
            (
                "Minv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Minv", dimMoles/dimMass, 0.0)
        )
    );

    volScalarField& Minv = tMinv();

    forAll(species_, j)
    {
        Minv += (Y_[j] + Z_[j])/dimensionedScalar("M", dimMass/dimMoles, M_[j]);
    }

    forAll(species_, j)
    {
        W_[j] == Y_[j]/Minv/dimensionedScalar("M", dimMass/dimMoles, M_[j]);
        X_[j] == Z_[j]/Minv/dimensionedScalar("M", dimMass/dimMoles, M_[j]);
    }
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::Ytot()
{
    tmp<volScalarField> tYtot
    (
        new volScalarField
        (
            IOobject
            (
                "Ytot",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Ytot", dimless, 0.0)
        )
    );

    volScalarField& Ytot = tYtot();

    forAll(species_, i)
    {
        Ytot += Y_[i];
    }

    return tYtot;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::Ztot()
{
    tmp<volScalarField> tZtot
    (
        new volScalarField
        (
            IOobject
            (
                "Ztot",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Ztot", dimless, 0.0)
        )
    );

    volScalarField& Ztot = tZtot();

    forAll(species_, i)
    {
        Ztot += Z_[i];
    }

    return tZtot;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::mVapor()
{
    tmp<volScalarField> tm
    (
        new volScalarField
        (
            IOobject
            (
                "m",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("m", dimMass, 0.0)
        )
    );

    tmp<volScalarField> tsumYm
    (
        new volScalarField
        (
            IOobject
            (
                "sumYm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("sumYm", dimless/dimMass, 0.0)
        )
    );

    volScalarField& sumYm = tsumYm();
    volScalarField& m = tm();

    tmp<volScalarField> tYtot = Ytot();
    volScalarField& Ytot = tYtot();

    Ytot.max(0.0);

    forAll(species_, j)
    {
        dimensionedScalar mj("m", dimMass, M()[j]*1E-03/6.022140857E23);

        sumYm += (Y_[j] / mj);
    }

    m == Ytot / sumYm;

    return tm;
}

Foam::scalar Foam::fluidThermo::getDiffusivity
(
    const label a,
    const label b,
    const scalar T,
    const scalar p
)
{
    if(diffusivityModels_.empty())
    {
        getDiffusivityModels();
    }

    label k = diffIndex_[a][b];

    return (k == -1) ? 0.0 : diffusivityModels_[k].value(T, p);
}

Foam::tmp<Foam::scalarField> Foam::fluidThermo::getDiffusivity
(
    const label a,
    const label b,
    const scalarField& T,
    const scalarField& p
)
{
    if(diffusivityModels_.empty())
    {
        getDiffusivityModels();
    }

    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& D = tD();

    forAll(T, i)
    {
        D[i] = this->getDiffusivity(a, b, T[i], p[i]);
    }

    return tD;
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::getDiffusivity
(
    const label a,
    const label b,
    const volScalarField& T,
    const volScalarField& p
)
{
    if(diffusivityModels_.empty())
    {
        getDiffusivityModels();
    }

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("D", dimLength*dimLength/dimTime, 0.0)
        )
    );

    volScalarField& D = tD();

    scalarField& Dcells = D.internalField();
    const scalarField& Tcells = T.internalField();
    const scalarField& pcells = p.internalField();

    forAll(Dcells, i)
    {
        Dcells[i] = this->getDiffusivity(a, b, Tcells[i], pcells[i]);
    }

    forAll(D.boundaryField(), patchi)
    {
        scalarField& Dpatch = D.boundaryField()[patchi];
        const scalarField& Tpatch = T.boundaryField()[patchi];
        const scalarField& ppatch = p.boundaryField()[patchi];

        forAll(Dpatch, i)
        {
            Dpatch[i] = this->getDiffusivity(a, b, Tpatch[i], ppatch[i]);
        }
    }

    return tD;
}

namespace Foam // Necessary for Doxygen
{

void fluidThermo::limitVaporDensity(surfaceScalarField& rho)
{
    rho.max(rhovMin_);
    rho.min(rhovMax_);
}

void fluidThermo::limitLiquidDensity(surfaceScalarField& rho)
{
    rho.max(rholMin_);
    rho.min(rholMax_);
}

void fluidThermo::limitVaporDensity(volScalarField& rho)
{
    rho.max(rhovMin_);
    rho.min(rhovMax_);
}

void fluidThermo::limitLiquidDensity(volScalarField& rho)
{
    rho.max(rholMin_);
    rho.min(rholMax_);
}

void fluidThermo::limitVaporDensity(scalarField& rho)
{
    rho = max(rho, rhovMin_);
    rho = min(rho, rhovMax_);
}

void fluidThermo::limitLiquidDensity(scalarField& rho)
{
    rho = max(rho, rholMin_);
    rho = min(rho, rholMax_);
}

} // End Foam namespace

bool Foam::fluidThermo::propertyFound
(
    const label i,
    const word propertyName,
    const phaseType phase
)
{
    HashTable< List<Switch> > *propertiesAvailablePtr(NULL);

    if (phase == VAPOR)
    {
        propertiesAvailablePtr = &propertiesAvailableVapor_;
    }
    else if (phase == LIQUID)
    {
        propertiesAvailablePtr = &propertiesAvailableLiquid_;
    }
    else
    {
        FatalErrorIn
        (
         "Foam::fluidThermo::propertyFound(const label i, const word propertyName, const phaseType phase)"
        )   << "No valid phase specified. Must be either vapor or liquid." << nl
            << exit(FatalError);
    }

    if (!propertiesAvailablePtr->found(propertyName))
    {
        return false;
    }
    else
    {
        List<Switch> propertiesAvailableList
            = propertiesAvailablePtr->find(propertyName)();

        if (!propertiesAvailableList[i])
        {
            return false;
        }
    }

    return true;
}

// ************************************************************************* //
