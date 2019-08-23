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

\file twoMomentLogNormal.H
\brief two-moment log-normal aerosol model class

Solves the PBE by assuming a log-normal distribution of which the width is
fixed. The distribution is closed by solving the number concentration transport
equation. Explicit right-hand side source terms are provided for the \f$Y_j\f$
and \f$Z_j\f$-equations.

*/

#ifndef twoMomentLogNormal_H
#define twoMomentLogNormal_H

#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                        Class twoMomentLogNormal Declaration
\*---------------------------------------------------------------------------*/

class twoMomentLogNormal
:
    public aerosolModel
{
    //- Private Member Functions

        //- Disallow copy construct
        twoMomentLogNormal(const twoMomentLogNormal&);

        //- Disallow default bitwise assignment
        void operator=(const twoMomentLogNormal&);


protected:

    //- Protected Data

        //- Number concentration field
        volScalarField M_;

        //- Number concentration rate change field
        volScalarField J_;

        //- Coalescence coefficient
        volScalarField f_;

        //- Mass concentration rate change field
        PtrList<volScalarField> I_;

        //- Geometric standard deviation
        scalar sigma_;


    //- Protected Member Functions

        //- Update the source terms
        void updateSources();

        //- Update the drift fields
        void updateDrift();


public:

    //- Runtime type information
    TypeName("twoMomentLogNormal");


    // Constructors

        //- Construct from components
        twoMomentLogNormal
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties
        );


    //- Destructor
    virtual ~twoMomentLogNormal();


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

} // End namespace twoMomentLogNormals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
