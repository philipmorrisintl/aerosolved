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

\file FullerShettlerGiddings.H
\brief Provides the Fuller-Schettler-Giddings diffusivity

Computes the binary diffusivity between a pair of vapor species based on the
Fuller-Schettler-Giddings relation, see \cite Fuller:1966aa. This relation
requires the 'diffusion volumes' to be specified for both species.

*/

#ifndef FullerShettlerGiddings_H
#define FullerShettlerGiddings_H

#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

// Forward declaration of classes
class FullerShettlerGiddings;

/*---------------------------------------------------------------------------*\
                        Class FullerShettlerGiddings Declaration
\*---------------------------------------------------------------------------*/

class FullerShettlerGiddings
:
    public diffusivityModel
{
    // Private data

        //- Diffusion volumes
        scalar Vda_;
        scalar Vdb_;

        //- Mole weights
        scalar Ma_;
        scalar Mb_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const FullerShettlerGiddings&);


public:

    //- Runtime type information
    TypeName("FullerShettlerGiddings");


    // Constructors

        FullerShettlerGiddings
        (
            const word& entryName,
            const dictionary& dict,
            const dictionary species,
            const label a,
            const label b
        );

        //- Copy constructor
        FullerShettlerGiddings(const FullerShettlerGiddings& fsg);

        //- Construct and return a clone
        virtual tmp<diffusivityModel> clone() const
        {
            return tmp<diffusivityModel>(new FullerShettlerGiddings(*this));
        }


    //- Destructor
    virtual ~FullerShettlerGiddings();


    // Member Functions

        // Manipulation


        // Evaluation

            //- Return FullerShettlerGiddings value
            scalar value(const scalar T, const scalar p) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
