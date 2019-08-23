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

/**

\file subGridDepositionModel.H
\brief Sub-grid deposition model

This class contains the solution algorithm for the sub-grid deposition velocity
model. It can be used by a boundary condition, which provides the scaled carrier
velocity, scaled gravity and scaled cell-centered velocity, to determine if a
collision with the wall will occur, and if so, what the velocity of the particle
at the wall would be.

*/

#ifndef subGridDepositionModel_H
#define subGridDepositionModel_H

#include "scalar.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class subGridDepositionModel Declaration
\*---------------------------------------------------------------------------*/

class subGridDepositionModel
{
private:

    // Private Data

        //- Gas velocity
        scalar u_;

        //- Gravity
        scalar g_;

        //- Initial particle velocity
        scalar v0_;

        //- Particle radius
        scalar r_;

        //- Maximum number of Newton iterations
        label maxIter_;

        //- Solver tolerance
        scalar tolerance_;


    // Private member functions

        //- Disallow default bitwise copy construct
        subGridDepositionModel(const subGridDepositionModel&);

        //- Disallow default bitwise assignment
        void operator=(const subGridDepositionModel&);


public:

    // Constructors

        //- Construct from u, g, v0, r, maxIter and tolerance
        subGridDepositionModel
        (
            const scalar u,
            const scalar g,
            const scalar v0,
            const scalar r,
            const scalar maxIter = 99,
            const scalar tolerance = 1e-8
        );

    //- Destructor
    virtual ~subGridDepositionModel();


    // Member Functions

        //- Check if the particle will hit the wall
        inline bool collision() const;

        //- Particle velocity given a time
        inline scalar x(const scalar t) const;

        //- Particle position given a time
        inline scalar v(const scalar t) const;

        //- Compute interception time
        inline scalar t() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "subGridDepositionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
