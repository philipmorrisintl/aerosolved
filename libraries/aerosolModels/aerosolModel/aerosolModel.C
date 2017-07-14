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

#include "aerosolModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aerosolModel, 0);
    defineRunTimeSelectionTable(aerosolModel, dictionary);

    #ifndef DOXYGENIGNORE

    template<>
    const char* NamedEnum
    <
        aerosolModelType,2
    >::names[] =
    {
        "moment",
        "sectional"
    };

    template<>
    const char* NamedEnum
    <
        aerosolSizeDistributionType,5
    >::names[] =
    {
        "none",
        "linear",
        "logarithmic",
        "geometric",
        "list"
    };

    template<>
    const char* NamedEnum
    <
        aerosolSizePositionType,3
    >::names[] =
    {
        "NA",
        "center",
        "interface"
    };

    template<>
    const char* NamedEnum
    <
        distMethod,3
    >::names[] =
    {
        "twoMoment",
        "fourMoment",
        "hybrid"
    };

    #endif
}

const Foam::NamedEnum<aerosolModelType, 2>
    Foam::aerosolModel::modelTypeNames;
const Foam::NamedEnum<aerosolSizeDistributionType, 5>
    Foam::aerosolModel::sizeDistributionTypeNames;
const Foam::NamedEnum<aerosolSizePositionType, 3>
    Foam::aerosolModel::sizePositionTypeNames;
const Foam::NamedEnum<distMethod, 3>
    Foam::aerosolModel::distMethodNames;

const Foam::word Foam::aerosolModel::dictName("aerosolProperties");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModel::aerosolModel
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
    thermoPtr_(fluidThermo::New(mesh)),
    M_(0),
    V_(0),
    S_(0),
    J_(0),
    phid_(0),
    HvapS_
    (
        IOobject
        (
            "Hvap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),
    y_(0),
    x_(0),
    P_(0),
    q_(0.0),
    yMin_(0.0),
    yMax_(0.0),
    params_(subDict("aerosolModelParameters")),
    doDrift_(false),
    doCoa_(false),
    doNuc_(false),
    doCond_(false),
    doCorrSizeDist_(false),
    doMonitors_(false),
    scalarMonitorPtrs_(0),
    vectorMonitorPtrs_(0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModel::~aerosolModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::fluidThermo& Foam::aerosolModel::thermo()
{
    return thermoPtr_();
}

const Foam::fluidThermo& Foam::aerosolModel::thermo() const
{
    return thermoPtr_();
}

const Foam::fvMesh& Foam::aerosolModel::mesh() const
{
    return mesh_;
}

bool Foam::aerosolModel::read()
{
    thermo().readProperty("Hvap", fluidThermo::LIQUID, thermo().speciesPhaseChange());

    params_.lookup("doDrift") >> doDrift_;
    params_.lookup("doCoa") >> doCoa_;
    params_.lookup("doNuc") >> doNuc_;
    params_.lookup("doCond") >> doCond_;
    params_.lookup("doCorrSizeDist") >> doCorrSizeDist_;
    params_.lookup("doMonitors") >> doMonitors_;

    Info << "Aerosol model: Drift is switched " << (doDrift_ ? "on" : "off") << endl;
    Info << "               Coalescence is switched " << (doCoa_ ? "on" : "off") << endl;
    Info << "               Nucleation is switched " << (doNuc_ ? "on" : "off") << endl;
    Info << "               Condensation is switched " << (doCond_ ? "on" : "off") << endl;
    Info << "               Size distribution correction is switched " << (doCorrSizeDist_ ? "on" : "off") << endl;
    Info << "               Monitors are switched " << (doMonitors_ ? "on" : "off") << endl;

    setMonitors();

    return true;
}

void Foam::aerosolModel::setMonitors()
{
    if (doMonitors_)
    {
        scalarMonitorPtrs_.clearStorage();
        vectorMonitorPtrs_.clearStorage();

        scalarMonitorPtrs_.insert
        (
            "Jcoa",
            new volScalarField
            (
                IOobject
                (
                    "monitor.Jcoa",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("Jcoa", dimless/dimVolume/dimTime, 0.0)
            )
        );

        scalarMonitorPtrs_.insert
        (
            "Jnuc",
            new volScalarField
            (
                IOobject
                (
                    "monitor.Jnuc",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("Jnuc", dimless/dimVolume/dimTime, 0.0)
            )
        );

        scalarMonitorPtrs_.insert
        (
            "znuc",
            new volScalarField
            (
                IOobject
                (
                    "monitor.znuc",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("znuc", dimMass, 0.0)
            )
        );

        scalarMonitorPtrs_.insert
        (
            "dcm",
            new volScalarField
            (
                IOobject
                (
                    "monitor.dcm",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("dcm", dimLength, 0.0)
            )
        );

        scalarMonitorPtrs_.insert
        (
            "dmm",
            new volScalarField
            (
                IOobject
                (
                    "monitor.dmm",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("dmm", dimLength, 0.0)
            )
        );

        forAll(thermo().species(), j)
        {
            scalarMonitorPtrs_.insert
            (
                word("S." + Foam::name(j)),
                new volScalarField
                (
                    IOobject
                    (
                        word("monitor.S." + Foam::name(j)),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("S", dimMass/dimVolume/dimTime, 0.0)
                )
            );
        }
    }
}

void Foam::aerosolModel::clearMonitors()
{
    clearScalarMonitors();
    clearVectorMonitors();
}

void Foam::aerosolModel::clearScalarMonitors()
{
    if (doMonitors_)
    {
        forAll(thermo().species(), j)
        {
            word name("S." + Foam::name(j));
            *scalarMonitorPtrs_[name] *= 0.0;
        }

        *scalarMonitorPtrs_["Jcoa"] *= 0.0;
        *scalarMonitorPtrs_["Jnuc"] *= 0.0;
        *scalarMonitorPtrs_["znuc"] *= 0.0;
        *scalarMonitorPtrs_["dcm"] *= 0.0;
        *scalarMonitorPtrs_["dmm"] *= 0.0;
    }
}

void Foam::aerosolModel::clearVectorMonitors()
{
    if (doMonitors_)
    {
    }
}

void Foam::aerosolModel::updateHvapS()
{

    PtrList< DataEntry<scalar> >& dataEntriesListHvap_l =  thermo().getProperty("Hvap", fluidThermo::LIQUID);

    const dictionary& species = thermo().species();

    // Internal field

    const scalarField& TCells = thermo().T().internalField();

    scalarField& HvapSCells = HvapS_.internalField();

    HvapSCells = 0.0;

    forAll(species, i)
    {
        word speciesName = species.keys()[i];

        if (thermo().propertyFound(i, "Hvap", fluidThermo::LIQUID))
        {
            forAll(HvapSCells, cellj)
            {
                HvapSCells[cellj] += S_[i].internalField()[cellj] * dataEntriesListHvap_l[i].value(TCells[cellj]);
            }
        }
    }

    forAll(HvapS_.boundaryField(), patchi)
    {
        scalarField& T = thermo().T().boundaryField()[patchi];
        scalarField& HvapSPatch = HvapS_.boundaryField()[patchi];

        HvapSPatch = 0.0;

        forAll(species, i)
        {
            word speciesName = species.keys()[i];

            const scalarField& S = S_[i].boundaryField()[patchi];

            if (thermo().propertyFound(i, "Hvap", fluidThermo::LIQUID))
            {
                forAll(HvapSPatch, facej)
                {
                    HvapSPatch[facej] += S[facej] * dataEntriesListHvap_l[i].value(T[facej]);
                }
            }
        }
    }
}

void Foam::aerosolModel::limitWallFlux
(
    Foam::surfaceScalarField& phi
) const
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& p = mesh_.boundaryMesh()[patchI];

        if (p.type() == "wall")
        {
            scalarField& phiWall = phi.boundaryField()[patchI];

            phiWall = max(phiWall, 0.0);
        }
    }
}

void Foam::aerosolModel::updateSizeDistribution()
{
    sizeDistType_ =
        sizeDistributionTypeNames.read(params_.lookup("sizeDistributionType"));

    switch(sizeDistType_)
    {
        case NOSIZEDIST:
        {
            yMin_ = 0.0;
            yMax_ = 0.0;
            P_ = 0;
            q_ = 0.0;

            y_.setSize(0);
            x_.setSize(0);
        }
        break;

        case LINEARSIZEDIST:
        {
            params_.lookup("yMin") >> yMin_;
            params_.lookup("yMax") >> yMax_;
            params_.lookup("P") >> P_;
            q_ = 0.0;

            if(yMin_ >= yMax_)
            {
                FatalErrorIn("Foam::aerosolModel::read()")
                    << "yMax must be larger than yMin." << exit(FatalError);
            }

            y_.setSize(P_+1);
            x_.setSize(P_);

            forAll(x_, i)
            {
                y_[i] = yMin_ + scalar(i)*(yMax_-yMin_)/P_;
                x_[i] = y_[i] + 0.5*(yMax_-yMin_)/P_;
            }

            y_[P_] = yMax_;
        }
        break;

        case LOGARITHMICSIZEDIST:
        {
            params_.lookup("yMin") >> yMin_;
            params_.lookup("yMax") >> yMax_;
            params_.lookup("P") >> P_;

            if(yMin_ >= yMax_ || yMin_ < 0.0 || yMax_ < 0.0)
            {
                FatalErrorIn("Foam::aerosolModel::read()")
                    << "yMin and yMax must be positive, and yMax must be larger than yMin." << exit(FatalError);
            }

            if (yMin_ == 0.0)
            {
                FatalErrorIn("Foam::aerosolModel::read()")
                    << "In a logarithmic distribution yMin cannot be zero." << exit(FatalError);
            }

            y_.setSize(P_+1);
            x_.setSize(P_);

            scalar a = pow(yMax_/yMin_, 1.0/P_);

            y_[0] = yMin_;
            x_[0] = yMin_ * pow(a, 0.5);

            for (label i = 1; i < P_; i++)
            {
                y_[i] = y_[i-1] * a;
                x_[i] = x_[i-1] * a;
            }

            y_[P_] = yMax_;
        }
        break;

        case GEOMETRICSIZEDIST:
        {
            params_.lookup("yMin") >> yMin_;
            params_.lookup("yMax") >> yMax_;
            params_.lookup("P") >> P_;
            params_.lookup("q") >> q_;

            if(yMin_ >= yMax_ || yMin_ < 0.0 || yMax_ < 0.0)
            {
                FatalErrorIn("Foam::aerosolModel::read()")
                    << "yMin and yMax must be positive, and yMax must be larger than yMin." << exit(FatalError);
            }

            if(q_ == 1.0)
            {
                FatalErrorIn("Foam::aerosolModel::read()")
                    << "q (the ratio between last and first cell) should be unequal to unity. Use a linear distribution for this." << exit(FatalError);
            }

            y_.setSize(P_+1);
            x_.setSize(P_);

            scalar r = pow(q_, 1.0/(scalar(P_)-1.0));

            scalar dy = (yMax_ - yMin_)*(1.0-r) / (1.0-q_*r);

            y_[0] = yMin_;
            y_[1] = yMin_ + dy;

            for (label i = 2; i < P_; i++)
            {
                y_[i] = y_[i-1]*(1+r) - r*y_[i-2];
            }

            y_[P_] = yMax_;

            forAll(x_, i)
            {
                x_[i] = (y_[i+1]+r*y_[i])/(r+1.0);
            }
        }
        break;

        case LIST:
        {
            params_.lookup("x") >> x_;
            params_.lookup("y") >> y_;

            if (x_.size() < 1)
            {
                FatalErrorIn("void Foam::aerosolModel::updateSizeDistribution()")
                    << "The size of x must be at least one." << exit(FatalError);
            }

            if ((x_.size()+1) != y_.size())
            {
                FatalErrorIn("void Foam::aerosolModel::updateSizeDistribution()")
                    << "The size of x must be the size of y minus one." << exit(FatalError);
            }

            for (label i = 0; i < (x_.size()-1); i++)
            {
                if (x_[i] > x_[i+1])
                {
                    FatalErrorIn("void Foam::aerosolModel::updateSizeDistribution()")
                        << "List x must be ordered from small to large." << exit(FatalError);
                }
            }

            for (label i = 0; i < (y_.size()-1); i++)
            {
                if (y_[i] > y_[i+1])
                {
                    FatalErrorIn("void Foam::aerosolModel::updateSizeDistribution()")
                        << "List y must be ordered from small to large." << exit(FatalError);
                }
            }

            P_ = x_.size();
        }
        break;

    }

    return;
}

// ************************************************************************* //
