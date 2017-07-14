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

\file compressible.H
\brief The fluidThermo compressible model

This is the main mixture model as discussed in \cite thesis. Following Eq. (2.8)
from that work, the density is written as \f$\rho=\psi p\f$ with \f$\psi\f$ a
function of composition, species-specific ideal gas laws for vapors and
species-specific densities of pure incompressible liquids. The updateRho member
function uses Eq. (2.8) \cite thesis to compute \f$\rho\f$. The updatePsi member
function uses the definition of \f$\psi\f$ to update \f$\psi\f$. Eq. (2.8) can
be used in the continuity equation to find an evolution equation for the
pressure, e.g., see Eq. (3.42) in \cite thesis, consisting of both a Laplacian
and partial time deriviative of pressure.

*/

#ifndef compressible_H
#define compressible_H

#include "fluidThermo.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermos
{

/*---------------------------------------------------------------------------*\
                         Class compressible Declaration
\*---------------------------------------------------------------------------*/

class compressible
:
    public fluidThermo
{
public:
    // Data types

        //- Cp or Cv identifier
        enum heatCapacityType{CP, CV};

private:
    // Private data


    // Private Member Functions

        //- Disallow copy construct
        compressible(const compressible&);

        //- Disallow default bitwise assignment
        void operator=(const compressible&);


protected:

    // Protected data

        const fvMesh& mesh_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("compressible");


    // Constructors

        //- Construct from components
        compressible(const fvMesh& mesh);


    //- Destructor
    virtual ~compressible();


    // Member Functions

        // Thermodynamics

            //- Update density, based on psi and pressure
            void updateRho();

            //- Update compressibility, based on the available pressure and
            //- density
            void updatePsi();

            //- Vapor density
            tmp<volScalarField> rhoY();

            //- Liquid density
            tmp<volScalarField> rhoZ();

            //- Vapor compressibility
            tmp<volScalarField> psiY();

            //- Current mixture density
            tmp<volScalarField> rhoMix();

            //- Update mixture heat capacity at constant volume
            void updateCvEff();

            //- Update mixture heat capacity at constant pressure
            void updateCpEff();

            //- Update given heat capacity field based on given heat capacity type
            void updateCxEff(volScalarField& CxEff, const heatCapacityType CpCv);

            //- Current liquid density
            tmp<volScalarField> rhoLiquid();

            //- Current vapor density
            tmp<volScalarField> rhoVapor();

        // I-O

            //- Read info from dictionaries
            bool read();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidThermos
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
