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

\file StokesEinstein.H
\brief The Stokes-Einstein droplet diffusivity model

See Eq. (7.8) and Sec. 7.2.2 in \cite thesis.

*/

#ifndef StokesEinstein_H
#define StokesEinstein_H

#include "brownianDiffusionModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace brownianDiffusionModels
{

/*---------------------------------------------------------------------------*\
                       Class StokesEinstein Declaration
\*---------------------------------------------------------------------------*/

class StokesEinstein
:
    public brownianDiffusionModel
{
    // Private data

        //- Critical diameter
        scalar Dcrit_;

    // Private Member Functions

        //- Disallow copy construct
        StokesEinstein(const StokesEinstein&);

        //- Disallow default bitwise assignment
        void operator=(const StokesEinstein&);


protected:


public:

    //- Runtime type information
    TypeName("StokesEinstein");


    // Constructors

        //- Construct from components
        StokesEinstein
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        );


    //- Destructor
    virtual ~StokesEinstein();


    // Member Functions

        // Thermodynamics

            //- Update diffusivity
            void update();

        // I-O

            //- Read info from dictionaries
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace brownianDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
