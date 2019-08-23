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

\file sectionalInterpolation.H
\brief Base class of the sectional interpolation scheme

This class determines how to interpolate sectional number concentrations between
representative sizes. It contains member functions which can add particles to
the distribution, provided a size of those particles.

*/

#ifndef sectionalInterpolation_H
#define sectionalInterpolation_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "sectionalDistribution.H"
#include "secIntData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sectionalInterpolation Declaration
\*---------------------------------------------------------------------------*/

class sectionalInterpolation
:
    public aerosolSubModelBase
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sectionalInterpolation(const sectionalInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const sectionalInterpolation&);


protected:

    // Protected Data

        //- Reference to the sectional distribution
        sectionalDistribution& distribution_;


public:

    //- Runtime type information
    TypeName("sectionalInterpolation");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             sectionalInterpolation,
             dictionary,
             (
                aerosolModel& aerosol,
                sectionalDistribution& distribution,
                const dictionary& dict
             ),
             (aerosol, distribution, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        sectionalInterpolation
        (
            const word& modelType,
            aerosolModel& aerosol,
            sectionalDistribution& distribution,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected sectional distribution
        static autoPtr<sectionalInterpolation> New
        (
            aerosolModel& aerosol,
            sectionalDistribution& distribution,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sectionalInterpolation();

    //- Access


    // Member Functions

        //- Compute the indices and weights given a size
        virtual secIntData interp(const scalar& s) const = 0;

        //- Add to M using the interpolation scheme
        inline void addToM
        (
            const scalar& s,
            const scalar& M,
            const label& celli
        );

        //- Add to M using a provided interpolation scheme
        inline void addToM
        (
            const secIntData& idata,
            const scalar& s,
            const scalar& M,
            const label& celli
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sectionalInterpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
