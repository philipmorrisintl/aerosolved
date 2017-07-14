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

\file compressibleVreman.H
\brief The compressibe Vreman LES model

LES model based on the work of Bert Vreman, see \cite Vreman:1994aa.

*/

#ifndef compressibleVreman_H
#define compressibleVreman_H

#include "turbModel.H"
#include "DataEntry.H"
#include "fvcGrad.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbModels
{

/*---------------------------------------------------------------------------*\
                    Class compressibleVreman Declaration
\*---------------------------------------------------------------------------*/

class compressibleVreman
:
    public turbModel
{
    // Private data

        scalar cs_;

        scalar c1_;

        scalar c2_;

        scalar Prt_;


    // Private Member Functions

        //- Disallow copy construct
        compressibleVreman(const compressibleVreman&);

        //- Disallow default bitwise assignment
        void operator=(const compressibleVreman&);


protected:


public:

    //- Runtime type information
    TypeName("compressibleVreman");


    // Constructors

        //- Construct from components
        compressibleVreman
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Destructor
    virtual ~compressibleVreman();

    // Member Functions

        // Thermodynamics

            //- Update effective turbulence
            void update();

       tmp<volScalarField> K(volVectorField& U);

        // I-O

            //- Read info from dictionaries
            bool read();


        inline const scalar& Prt() const
        {
             return Prt_;
        }

        inline const scalar& c1() const
        {
             return c1_;
        }

        inline const scalar& c2() const
        {
             return c2_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
