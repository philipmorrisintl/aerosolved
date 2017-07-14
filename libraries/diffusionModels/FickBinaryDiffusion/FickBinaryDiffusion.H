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

/**

\file FickBinaryDiffusion.H
\brief Fick binary diffusion model

Simple model where the effective diffusivity \f$\mathbb{D}_j^v\f$ is taken equal
to the binary diffusivity of species a in species b. In a binary system this
leads to consistency as both species diffuse with an equal velocity, if the
gradients of both species are equal and opposite to each other. In the case of a
non-diffusing aerosol, this requirement may not necessarily be guaranteed since
gradients may differ due to the presence of a droplet phase.

*/

#ifndef FickBinaryDiffusion_H
#define FickBinaryDiffusion_H

#include "diffusionModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{

/*---------------------------------------------------------------------------*\
                    Class FickBinaryDiffusion Declaration
\*---------------------------------------------------------------------------*/

class FickBinaryDiffusion
:
    public diffusionModel
{

    // Private Member Functions

        //- Disallow copy construct
        FickBinaryDiffusion(const FickBinaryDiffusion&);

        //- Disallow default bitwise assignment
        void operator=(const FickBinaryDiffusion&);

protected:

    // Protected data

        //- Name of inert species
        word inert_;

        //- Index in species dictionary pointing to inert species
        label inertIndex_;

public:

    //- Runtime type information
    TypeName("FickBinaryDiffusion");


    // Constructors

        //- Construct from components
        FickBinaryDiffusion
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~FickBinaryDiffusion();


    // Member Functions

        // Thermodynamics

            //- Update effective diffusivity
            void update();

        // I-O

            //- Read info from dictionaries
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
