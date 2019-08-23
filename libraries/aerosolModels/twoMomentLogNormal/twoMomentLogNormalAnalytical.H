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

\file twoMomentLogNormalAnalytical.H
\brief Analytical extension of the two-moment log-normal aerosol model

A copy of the twoMomentLogNormal model, but does not provide explicit source
terms for the \f$Y_j\f$ and \f$Z_j\f$ equations because these source terms are
solved analytically in the solvePost() step. Generally, the
twoMomentLogNormalAnalytical moment is more stable than the twoMomentLogNormal
model, and is therefore recommended to be used

*/

#ifndef twoMomentLogNormalAnalytical_H
#define twoMomentLogNormalAnalytical_H

#include "aerosolModel.H"
#include "twoMomentLogNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
               Class twoMomentLogNormalAnalytical Declaration
\*---------------------------------------------------------------------------*/

class twoMomentLogNormalAnalytical
:
    public twoMomentLogNormal
{
    //- Private Member Functions

        //- Disallow copy construct
        twoMomentLogNormalAnalytical(const twoMomentLogNormalAnalytical&);

        //- Disallow default bitwise assignment
        void operator=(const twoMomentLogNormalAnalytical&);


public:

    //- Runtime type information
    TypeName("twoMomentLogNormalAnalytical");


    // Constructors

        //- Construct from components
        twoMomentLogNormalAnalytical
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties
        );


    //- Destructor
    virtual ~twoMomentLogNormalAnalytical();


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

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoMomentLogNormalAnalyticals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
