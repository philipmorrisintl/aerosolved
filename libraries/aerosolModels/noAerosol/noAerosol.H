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

#ifndef noAerosol_H
#define noAerosol_H

#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                        Class noAerosol Declaration
\*---------------------------------------------------------------------------*/

class noAerosol
:
    public aerosolModel
{
    //- Private Member Functions

        //- Disallow copy construct
        noAerosol(const noAerosol&);

        //- Disallow default bitwise assignment
        void operator=(const noAerosol&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noAerosol
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties
        );


    //- Destructor
    virtual ~noAerosol();


    // Member Functions

        //- Correct the aerosol model
        virtual void correctModel();

        //- Solution step before the mass fraction solution
        virtual void solvePre();

        //- Solution step after the mass fraction solution
        virtual void solvePost();

        //- Right-hand side source term
        virtual tmp<fvScalarMatrix> R(const volScalarField& Y) const;

        //- Heat release rate
        virtual tmp<volScalarField> Qdot() const;

        //- Moment p and q ratio mean diameter
        virtual tmp<volScalarField> meanDiameter
        (
            const scalar p,
            const scalar q
        ) const;

        //- Moment p median diameter
        virtual tmp<volScalarField> medianDiameter
        (
            const scalar p
        ) const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noAerosols
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
