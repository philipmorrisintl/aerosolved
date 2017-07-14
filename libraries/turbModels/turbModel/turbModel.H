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

\file turbModel.H
\brief Turbulence modeling

This library provides the sub-grid turbulence modeling. The five key functions
of this model are:

- update(): update the model

- divDevRhoBeff(): return the extra NS term due to the filtered NS equations

- K(): extra energy production in temperature equation

- muEff() and kEff(): effective momentum and heat diffusion coefficients

- muTurb() and kTurb(): turbulent momentum and heat diffusion coefficients

*/

#ifndef turbModel_H
#define turbModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvmSup.H"
#include "fluidThermo.H"
#include "IOdictionary.H"

#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       class turbModel Declaration
\*---------------------------------------------------------------------------*/

class turbModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to thermo object
        fluidThermo& thermo_;

        //- Coefficients
        dictionary coeffs_;

        //- Name of velocity field
        word velocityName_;

        //- Name of density field
        word densityName_;

        //- Name of dynamic mass flux field
        word massFluxName_;

        //- Name of laminar viscosity field
        word viscosityName_;

        //- Name of laminar diffusivity field
        word diffusivityName_;

        //- Effective turbulence viscosity field
        volScalarField muTurb_;

        //- Effective turbulence diffusivity field
        volScalarField kTurb_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to mass flux field
        const surfaceScalarField& phi_;

        //- Reference to laminar viscosity field
        const volScalarField& muEff_;

        //- Reference to laminar diffusivity field
        const volScalarField& kEff_;

        volScalarField delta_;

        scalar deltaCoeff_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        turbModel(const turbModel&);

        //- Disallow default bitwise assignment
        void operator=(const turbModel&);


public:

    //- Runtime type information
    TypeName("turbModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        turbModel,
        dictionary,
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        ),
        (mesh, thermo)
    );


    // Constructors

        //- Construct from components
        turbModel
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Selector
    static autoPtr<turbModel> New
    (
        const fvMesh& mesh,
        fluidThermo& thermo
    );


    //- Destructor
    virtual ~turbModel();


    // Member Functions

        // Access

            //- Turbulent viscosity
            inline const volScalarField& muTurb() const;

            //- Turbulent heat conductivity
            inline const volScalarField& kTurb() const;

            //- Effective viscosity
            inline const volScalarField& muEff() const;

            //- Effective heat conductivity
            inline const volScalarField& kEff() const;

            //- Return access to the thermo object
            inline fluidThermo& thermo();

            //- Update effective turbulence (pure virtual)
            virtual void update() = 0;

            //- Returns div(rho*dev(B)), the additional term due to the filtering of the NSE.
            tmp<fvVectorMatrix> divDevRhoBeff(volVectorField& U) const;

            //- K scalar field
            virtual tmp<volScalarField> K(volVectorField& U) = 0;


        // I-O

            //- Read info from dictionaries (pure virtual)
            virtual bool read() = 0;

            //- Return access to U
            inline const volVectorField& U() const;

            //- Return access to rho
            inline const volScalarField& rho() const;

            //- Return access to the delta field
            inline const volScalarField& delta() const;

            //- Return access to deltaCoeff scalar
            inline const scalar& deltaCoeff() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "turbModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
