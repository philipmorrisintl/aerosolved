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

#include "FickBinaryDiffusion.H"
#include "makeDiffusionTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
    makeDiffusionModelTypes(FickBinaryDiffusion, diffusionModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModels::FickBinaryDiffusion::FickBinaryDiffusion
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    diffusionModel(mesh, thermo),
    inert_(""),
    inertIndex_(-1)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionModels::FickBinaryDiffusion::~FickBinaryDiffusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionModels::FickBinaryDiffusion::update()
{
    volScalarField& DYinert = DY_[inertIndex_];

    const volScalarField& T = thermo_.T();
    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    DYinert *= 0;

    forAll(thermo_.species(), j)
    {
        if (j != inertIndex_)
        {
            DY_[j] == thermo_.getDiffusivity(j, inertIndex_, T, p1+p0);

            DYinert += DY_[j];
        }
    }

    // The inert species diffusivity is set to the average. This is a random
    // choice, but works for a binary system, to which Fick actually is limited.

    DYinert /= (thermo_.nSpecies()-1);
}

bool Foam::diffusionModels::FickBinaryDiffusion::read()
{
    coeffs_.lookup("inertSpecies") >> inert_;

    if (thermo_.species().found(inert_))
    {
        forAll(thermo_.species(), j)
        {
            if (thermo_.species().keys()[j] == inert_)
            {
                inertIndex_ = j;
                break;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::diffusionModels::FickBinaryDiffusion::read()"
        )   << "Could not find inert species in dictionary" << nl
            << exit(FatalError);
    }

    return true;
}


// ************************************************************************* //
