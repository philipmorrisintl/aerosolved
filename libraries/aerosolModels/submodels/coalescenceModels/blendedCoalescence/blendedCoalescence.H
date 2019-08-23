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

#ifndef blendedCoalescence_H
#define blendedCoalescence_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class blendedCoalescence Declaration
\*---------------------------------------------------------------------------*/

class blendedCoalescence
:
    public coalescenceModel
{
    // Private Data

        //- Low Knudsen model
        autoPtr<coalescenceModel> coaModel1_;

        //- High Knudsen model
        autoPtr<coalescenceModel> coaModel2_;


    // Private member functions

        //- Disallow default bitwise copy construct
        blendedCoalescence(const blendedCoalescence&);

        //- Disallow default bitwise assignment
        void operator=(const blendedCoalescence&);

        //- Blend two coaData objects given Knudsen
        coaData blend
        (
            const coaData& coa1,
            const coaData& coa2,
            const scalar phi
        ) const;


public:

    //- Runtime type information
    TypeName("blended");


    // Constructors

        //- Construct from aerosol model
        blendedCoalescence(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~blendedCoalescence();


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
