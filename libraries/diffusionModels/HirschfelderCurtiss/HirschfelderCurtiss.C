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

#include "HirschfelderCurtiss.H"
#include "makeDiffusionTypes.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
    makeDiffusionModelTypes(HirschfelderCurtiss, diffusionModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModels::HirschfelderCurtiss::HirschfelderCurtiss
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    diffusionModel(mesh, thermo)
{
   read();
   update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionModels::HirschfelderCurtiss::~HirschfelderCurtiss()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionModels::HirschfelderCurtiss::update()
{
    thermo_.updateMoleFractions();

    const volScalarField& T = thermo_.T();
    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();
    const volScalarField& rho = thermo_.rho();
    const PtrList<volScalarField>& W = thermo_.W();
    const PtrList<volScalarField>& Y = thermo_.Y();

    // Determine D-star

    tmp<volScalarField> tSumYk
    (
        new volScalarField
        (
            IOobject
            (
                "sumYk",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("sumYk", dimless, 0.0)
        )
    );

    tmp<volScalarField> tSumWkOverDjk
    (
        new volScalarField
        (
            IOobject
            (
                "sumWkOverDjk",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("sumWkOverDjk", dimTime/sqr(dimLength), 0.0)
        )
    );

    volScalarField& sumYk = tSumYk();
    volScalarField& sumWkOverDjk = tSumWkOverDjk();

    forAll(thermo_.species(), j)
    {
        sumYk *= 0;
        sumWkOverDjk *= 0;

        forAll(thermo_.species(), k)
        {
            if (j != k)
            {
                volScalarField Djk = thermo_.getDiffusivity(j, k, T, p1+p0);

                sumYk += max(Y[k], 1E-10);
                sumWkOverDjk += max(W[k], 1E-10)/Djk;
            }
        }

        DY_[j] = sumYk/sumWkOverDjk;
    }

    // Compute interpolated sum of mass fractions

    surfaceScalarField sumY
    (
        IOobject
        (
            "sumY",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("sumY", dimless, 0.0)
    );

    forAll(thermo_.species(), j)
    {
        sumY += fvc::interpolate(Y[j]);
    }

    // Determine correction velocity, with snGrad(Yj)

    phic_ *= 0;

    forAll(thermo_.species(), j)
    {
        phic_ += fvc::interpolate(rho*DY_[j]) * fvc::snGrad(Y[j]) * mesh_.magSf();
    }

    phic_ /= sumY;
}

bool Foam::diffusionModels::HirschfelderCurtiss::read()
{
    return true;
}

// ************************************************************************* //
