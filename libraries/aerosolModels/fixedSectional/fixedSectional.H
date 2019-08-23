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

\file fixedSectional.H
\brief Fixed sectional aerosol model

Solves the PBE using a sectional discretization, in which the sections
(specified in terms of particle mass) are fixed in time and space. The
fixedSectional object relies on the fixedSectionalSystem, which, in turn,
provides the sectional distribution and sectional interpolation functionality

*/

#ifndef fixedSectional_H
#define fixedSectional_H

#include "aerosolModel.H"
#include "fixedSectionalSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                        Class fixedSectional Declaration
\*---------------------------------------------------------------------------*/

class fixedSectional
:
    public aerosolModel
{
    //- Private Member Functions

        //- Disallow copy construct
        fixedSectional(const fixedSectional&);

        //- Disallow default bitwise assignment
        void operator=(const fixedSectional&);


protected:

    //- Protected Data

        //- The sectional data
        autoPtr<fixedSectionalSystem> system_;

        //- Number concentration rate change field (monitor)
        volScalarField J_;

        //- Mass concentration rate change field (monitor)
        PtrList<volScalarField> I_;


    //- Protected Member Functions

        //- Update the drift fields
        void updateDrift();

        //- Solve the spatial part of the sectional mass fraction equations
        void solveSpatial();

        //- Solve the internal part of the sectional mass fraction equations
        void solveInternal();


public:

    //- Runtime type information
    TypeName("fixedSectional");


    // Constructors

        //- Construct from components
        fixedSectional
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties
        );


    //- Destructor
    virtual ~fixedSectional();


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
        virtual tmp<volScalarField> medianDiameter(const scalar p) const;

        //- Clear the condensation and nucleation rates
        virtual void clearRates();

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fixedSectionals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
