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

\file HirschfelderCurtiss.H
\brief Hirschfelder-Curtiss multi-species diffusion treatment

The Hirschfelder-Curtiss (see \cite Curtiss:1949aa) offers a consistent
multi-species vapor diffusion framework, by computing an effective diffusivity
\f$\mathbb{D}_j^v\f$ and a corrective flux which guarantees that the sum of all
mass fractions remains unity. Also, see \cite Giovangigli:1991aa.

*/

#ifndef HirschfelderCurtiss_H
#define HirschfelderCurtiss_H

#include "diffusionModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{

/*---------------------------------------------------------------------------*\
                      Class HirschfelderCurtiss Declaration
\*---------------------------------------------------------------------------*/

class HirschfelderCurtiss
:
    public diffusionModel
{

    // Private Member Functions

        //- Disallow copy construct
        HirschfelderCurtiss(const HirschfelderCurtiss&);

        //- Disallow default bitwise assignment
        void operator=(const HirschfelderCurtiss&);


protected:


public:

    //- Runtime type information
    TypeName("HirschfelderCurtiss");


    // Constructors

        //- Construct from components
        HirschfelderCurtiss
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~HirschfelderCurtiss();


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
