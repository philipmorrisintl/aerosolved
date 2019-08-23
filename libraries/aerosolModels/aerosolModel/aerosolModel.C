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

#include "aerosolModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aerosolModel, 0);
    defineRunTimeSelectionTable(aerosolModel, dictionary);
}

const Foam::word Foam::aerosolModel::aerosolPropertiesName
(
    "aerosolProperties"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModel::aerosolModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    IOdictionary
    (
        IOobject
        (
            aerosolProperties,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    thermo_(mesh),
    turbulencePtr_(),
    mesh_(mesh),
    coeffs_(modelType == "none" ? *this : subDict(modelType + "Coeffs")),
    modelType_(modelType),
    outputPropertiesPtr_(),
    condensation_(),
    nucleation_(),
    coalescence_(),
    drift_(),
    dMin_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("min"))
    ),
    dMax_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("max"))
    ),
    phiDrift_
    (
        IOobject
        (
            "phiDrift",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    tauDrift_
    (
        IOobject
        (
            "tauDrift",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
            "tau",
            dimMass/sqr(dimTime)/dimLength,
            symmTensor::zero
        )
    ),
    phiInertial_
    (
        IOobject
        (
            "phiInertial",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    phiBrownian_
    (
        IOobject
        (
            "phiBrownian",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    DDisp_
    (
        IOobject
        (
            "DDisp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("D", dimArea/dimTime, 0.0)
    ),
    DCont_(thermo_.contSpecies().size()),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffs_.lookupOrDefault<scalar>("residualAlpha", 1E-12)
    )
{
    read();

    if (!outputPropertiesPtr_.valid())
    {
        const fileName uniformPath(word("uniform")/"aerosolModels");

        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "outputProperties",
                    mesh_.time().timeName(),
                    uniformPath,
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    if (modelType != "none")
    {
        condensation_ =
            condensationModel::New
            (
                *this,
                subDict("submodels").subDict("condensation")
            );

        nucleation_ =
            nucleationModel::New
            (
                *this,
                subDict("submodels").subDict("nucleation")
            );

        coalescence_ =
            coalescenceModel::New
            (
                *this,
                subDict("submodels").subDict("coalescence")
            );
    }

    drift_.set(
        new driftFluxModel(*this, subDict("submodels"))
    );

    forAll(thermo_.contSpecies(), j)
    {
        DCont_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("DCont", thermo_.contSpecies()[j]),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DCont", dimArea/dimTime, 0.0)
            )
        );
    }

    forAll(thermo_.contSpecies(), j)
    {
        fields_.add(thermo_.Y()[j]);
    }

    forAll(thermo_.dispSpecies(), j)
    {
        fields_.add(thermo_.Z()[j]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModel::~aerosolModel()
{
    if (turbulencePtr_)
    {
        turbulencePtr_ = 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModel::correct()
{
    correctModel();

    if (this->drift().diffusion().type() != "none")
    {
        forAll(thermo_.contSpecies(), j)
        {
            DCont_[j] = this->drift().diffusion().D(j);
        }
    }
    else
    {
        forAll(thermo_.contSpecies(), j)
        {
            DCont_[j] *= 0.0;
        }
    }

    if
    (
        this->drift().diffusion().type() != "none"
     || this->drift().Brownian().type() != "none"
     || this->drift().inertial().type() != "none"
    )
    {
        phiDrift_ = this->drift().phi
        (
            phiInertial_,
            phiBrownian_,
            DDisp_,
            DCont_
        );

        tauDrift_ = this->drift().tau
        (
            phiInertial_,
            phiBrownian_,
            DDisp_,
            DCont_
        );
    }

    mvPhi_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            mesh_,
            fields_,
            this->phi(),
            mesh_.divScheme("div(mvConv)")
        );

    mvPhiDrift_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            mesh_,
            fields_,
            phiDrift_,
            mesh_.divScheme("div(mvConv)")
        );

    mvPhiInertial_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            mesh_,
            fields_,
            phiInertial_,
            mesh_.divScheme("div(mvConv)")
        );

    mvPhiBrownian_ =
        fv::convectionScheme<Foam::scalar>::New
        (
            mesh_,
            fields_,
            phiBrownian_,
            mesh_.divScheme("div(mvConv)")
        );
}

bool Foam::aerosolModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}

Foam::tmp<Foam::scalarField> Foam::aerosolModel::getRDeltaT()
{
    if(mesh_.foundObject<volScalarField>("rDeltaT"))
    {
        return tmp<scalarField>
        (
            mesh_.lookupObject<volScalarField>("rDeltaT")
        );
    }
    else
    {
        const scalarField& rho = this->rho().field();

        if(!rDeltaT_.valid())
        {
            rDeltaT_.set(new scalarField(rho.size()));
        }

        rDeltaT_() = 1.0/mesh_.time().deltaTValue();

        return tmp<scalarField>(rDeltaT_());
    }
}

// ************************************************************************* //
