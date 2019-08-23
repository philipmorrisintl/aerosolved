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

#include "Vreman.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Vreman<BasicTurbulenceModel>::correctNut()
{
    const volTensorField alpha(fvc::grad(this->U_));

    volTensorField beta
    (
        IOobject
        (
            "beta",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("beta", sqr(dimVelocity), tensor::zero)
    );

    beta = Foam::sqr(this->delta()) * (alpha.T() & alpha);

    const volScalarField B
    (
        beta.component(tensor::XX)*beta.component(tensor::YY)
      - beta.component(tensor::XY)*beta.component(tensor::XY)
      + beta.component(tensor::XX)*beta.component(tensor::ZZ)
      - beta.component(tensor::XZ)*beta.component(tensor::XZ)
      + beta.component(tensor::YY)*beta.component(tensor::ZZ)
      - beta.component(tensor::YZ)*beta.component(tensor::YZ)
    );

    const dimensionedScalar alphaSqrSmall
    (
        "small",
        alpha.dimensions()*alpha.dimensions(),
        SMALL
    );

    const dimensionedScalar BSmall
    (
        "small",
        B.dimensions(),
        SMALL
    );

    this->nut_ =
        2.5*sqr(Cs_)
      * Foam::sqrt(max(B,BSmall)/max(alpha && alpha, alphaSqrSmall));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Vreman<BasicTurbulenceModel>::Vreman
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    Smagorinsky<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.094
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Vreman<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void Vreman<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
