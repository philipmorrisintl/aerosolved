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

#ifndef coupledNucleation_H
#define coupledNucleation_H

#include "nucleationModel.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupledNucleation Declaration
\*---------------------------------------------------------------------------*/

class coupledNucleation
:
    public nucleationModel
{
    // Private Data

        //- Inner and outer iteration tolerance
        scalar TOL_;

        //- Max number of inner and outer iterations
        label maxInnerIter_;
        label maxOuterIter_;

    // Private member functions

        //- Compute the critical cluster composition
        void computeComposition
        (
            scalarList& w,
            scalarList& pm,
            scalar& alpha,
            const scalar& f,
            const scalar& Ke,
            const scalar& p,
            const scalarList& pVap,
            const scalarList& pSat,
            const scalarList& D,
            const scalarList& v
        ) const;

        //- Disallow default bitwise copy construct
        coupledNucleation(const coupledNucleation&);

        //- Disallow default bitwise assignment
        void operator=(const coupledNucleation&);


public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from aerosol model
        coupledNucleation(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~coupledNucleation();


    // Member Functions

        // Evolution

            //- Compute the nucleation data
            virtual nucData rate
            (
                const scalar& p,
                const scalar& T,
                const scalarList& Y,
                const scalarList& pSat,
                const scalarList& D,
                const scalarList& rhoDisp,
                const scalarList& sigma
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
