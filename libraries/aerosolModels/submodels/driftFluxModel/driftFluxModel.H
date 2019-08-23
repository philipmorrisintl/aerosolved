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

\file driftFluxModel.H
\brief Base class of the drift flux modeling

This class contains the diffusion, Brownian motion and inertial drift objects.
It also implements the computation of the drift flux stress tensor and
corrective velocity using the tau() and V() member functions, respectively.

*/

#ifndef driftFluxModel_H
#define driftFluxModel_H

#include "aerosolSubModelBase.H"
#include "diffusionModel.H"
#include "BrownianModel.H"
#include "inertialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class driftFluxModel Declaration
\*---------------------------------------------------------------------------*/

class driftFluxModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Data

        //- Diffusion model
        autoPtr<diffusionModel> diffusion_;

        //- Brownian model
        autoPtr<BrownianModel> Brownian_;

        //- Inertial model
        autoPtr<inertialModel> inertial_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        driftFluxModel(const driftFluxModel&);

        //- Disallow default bitwise assignment
        void operator=(const driftFluxModel&);


public:

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        driftFluxModel
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~driftFluxModel();


    // Access

        //- Sccess to the diffusion model
        virtual inline diffusionModel& diffusion()
        {
            return diffusion_();
        }

        //- Const access to the diffusion model
        virtual inline const diffusionModel& diffusion() const
        {
            return diffusion_();
        }

        //- Access to the brownian model
        virtual inline BrownianModel& Brownian()
        {
            return Brownian_();
        }

        //- Const access to the brownian model
        virtual inline const BrownianModel& Brownian() const
        {
            return Brownian_();
        }

        //- Access to the inertial model
        virtual inline inertialModel& inertial()
        {
            return inertial_();
        }

        //- Const access to the inertial model
        virtual inline const inertialModel& inertial() const
        {
            return inertial_();
        }


    // Member Functions

        //- Compute and return a diffusive flux given a mass fraction and
        //  diffusivity
        virtual tmp<surfaceScalarField> phiDiffusion
        (
            const volScalarField& Y,
            const volScalarField& D
        ) const;

        //- Compute and return a diffusive velocity given a mass fraction and
        //  diffusivity
        virtual tmp<volVectorField> VDiffusion
        (
            const volScalarField& Y,
            const volScalarField& D
        ) const;

        //- Compute and return the corrective velocity
        virtual tmp<volVectorField> V
        (
            const volVectorField& VInertial,
            const volVectorField& VBrownian,
            const volScalarField& DDisp,
            const PtrList<volScalarField>& DCont
        ) const;

        //- Compute and return the corrective flux
        virtual tmp<surfaceScalarField> phi
        (
            const surfaceScalarField& phiInertial,
            const surfaceScalarField& phiBrownian,
            const volScalarField& DDisp,
            const PtrList<volScalarField>& DCont
        ) const;

        //- Compute and return the mixture drift stress tensor
        virtual tmp<volSymmTensorField> tau
        (
            const surfaceScalarField& phiInertial,
            const surfaceScalarField& phiBrownian,
            const volScalarField& DDisp,
            const PtrList<volScalarField>& DCont
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
