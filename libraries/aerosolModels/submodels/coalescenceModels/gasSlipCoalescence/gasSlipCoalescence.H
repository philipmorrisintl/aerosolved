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

#ifndef gasSlipCoalescence_H
#define gasSlipCoalescence_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gasSlipCoalescence Declaration
\*---------------------------------------------------------------------------*/

class gasSlipCoalescence
:
    public coalescenceModel
{
    // Private Data

        //- A constant
        scalar A_;

        //- p coefficients
        scalarList p_;

        //- q coefficients
        scalarList q_;


    // Private member functions

        //- Disallow default bitwise copy construct
        gasSlipCoalescence(const gasSlipCoalescence&);

        //- Disallow default bitwise assignment
        void operator=(const gasSlipCoalescence&);


public:

    //- Runtime type information
    TypeName("gasSlip");


    // Constructors

        //- Construct from aerosol model
        gasSlipCoalescence(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~gasSlipCoalescence();


    // Member Functions

        // Evolution

            //- Compute the coalescence data
            virtual coaData rate
            (
                const scalar& p,
                const scalar& T,
                const scalar& mu,
                const scalar& rhog,
                const scalar& rhol,
                const scalar& d
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
