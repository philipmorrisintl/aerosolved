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

\file incompressible.H
\brief The fluidThermo incompressible model

In the incompressible model the density will stay constant in time. The density
is initially computed from mass fractions and their corresponding densities of
pure consistituent. The density is kept constant thereafter. This means that the
partial time derivative \f$\partial_t \rho=0\f$ inside the continuity equation, and
can therefore be left out of the pressure equation.

*/

#ifndef incompressible_H
#define incompressible_H

#include "fluidThermo.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermos
{

/*---------------------------------------------------------------------------*\
                         Class incompressible Declaration
\*---------------------------------------------------------------------------*/

class incompressible
:
    public fluidThermo
{
    // Private data

        //- Fixed uniform density value
        scalar rhoValue_;


    // Private Member Functions

        //- Disallow copy construct
        incompressible(const incompressible&);

        //- Disallow default bitwise assignment
        void operator=(const incompressible&);


protected:

    // Protected data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("incompressible");


    // Constructors

        //- Construct from components
        incompressible(const fvMesh& mesh);


    //- Destructor
    virtual ~incompressible();


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

            //- Update mixture heat capacity at constant volume, based on the
            //- available temperature and species vapor/liquid mass fractions,
            //- and the provided temperature-dependent definitions of the
            //- heat capacity at constant volume
            void updateCvEff();

            //- Update mixture heat capacity at constant pressure. The CvEff
            //- field is leading. CpEff is computed with CpEff=gamma*CvEff
            void updateCpEff();

            //- Mixture heat capacity for a patch, based on given temperature
            tmp<scalarField> CvEff
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Mixture heat capacity for a cell-set, based on given
            //- temperature. If no cell-set is given, internalField() is
            //- assumed
            tmp<scalarField> CvEff
            (
                const scalarField& T,
                labelList cells = labelList(0)
            ) const;

            //- Mixture heat capacity for a volScalarField, based on given temperature
            tmp<volScalarField> CvEff(const volScalarField& T) const;

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

#include "incompressibleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
