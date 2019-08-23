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

\file aerosolModel.H
\brief Main aerosol model base class

This class contains the implementation of different aerosol models, which are
responsible for the modeling of the Population Balance Equation (PBE).
Ultimately, the goal of solving the PBE (either in its true form, or in terms of
its moments) is to produce size information on the distribution of particles,
locally. In turn, this size information provides sub-models such as the
condensationModel and the driftFluxModel with relevant average diameters. The
aerosolModel class is the most important class in AeroSolved, and also contains
the aerosolThermo class, which describes the thermodynamical state of the
two-phase (i.e., continuous phase and dispersed phase) mixture.

*/

#ifndef aerosolModel_H
#define aerosolModel_H

#include "autoPtr.H"
#include "aerosolThermo.H"
#include "condensationModel.H"
#include "nucleationModel.H"
#include "coalescenceModel.H"
#include "driftFluxModel.H"
#include "IOdictionary.H"
#include "turbulentFluidThermoModel.H"
#include "multivariateScheme.H"
#include "convectionScheme.H"
#include "ListFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class aerosolModel Declaration
\*---------------------------------------------------------------------------*/

class aerosolModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow copy construct
        aerosolModel(const aerosolModel&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolModel&);


        //- the local time-scale
        autoPtr<scalarField> rDeltaT_;

protected:

    // Protected data

        //- Thermo model
        aerosolThermo thermo_;

        //- Reference to the turbulence model
        compressibleTurbulenceModel* turbulencePtr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;

        //- Condensation model
        autoPtr<condensationModel> condensation_;

        //- Nucleation model
        autoPtr<nucleationModel> nucleation_;

        //- Coalescence model
        autoPtr<coalescenceModel> coalescence_;

        //- Drift flux model
        autoPtr<driftFluxModel> drift_;

        //- Multivariate convection scheme for phi
        tmp<fv::convectionScheme<scalar>> mvPhi_;

        //- Multivariate convection scheme for phiDrift
        tmp<fv::convectionScheme<scalar>> mvPhiDrift_;

        //- Multivariate convection scheme for phiInertial
        tmp<fv::convectionScheme<scalar>> mvPhiInertial_;

        //- Multivariate convection scheme for phiBrownian
        tmp<fv::convectionScheme<scalar>> mvPhiBrownian_;

        //- Table of fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

        //- Minimum allowable diameter
        scalar dMin_;

        //- Maximum allowable diameter
        scalar dMax_;

        //- Corrective mass fraction drift flux [kg/s]
        surfaceScalarField phiDrift_;

        //- Mixture drift stress tensor [kg/m/s^2]
        volSymmTensorField tauDrift_;

        //- Inertial dispersed mass fraction drift flux [kg/s]
        surfaceScalarField phiInertial_;

        //- Brownian dispersed mass fraction drift flux [kg/s]
        surfaceScalarField phiBrownian_;

        //- Dispersed mass fraction diffusivity [m^2/s]
        volScalarField DDisp_;

        //- Continuous mass fraction diffusivities [m^2/s]
        PtrList<volScalarField> DCont_;

        //- Alpha residual for numerical stability
        dimensionedScalar residualAlpha_;

        //- List of effective particle number flux fields, used as monitor
        PtrList<surfaceScalarField> phiEff_;


public:

    //- Runtime type information
    TypeName("aerosolModel");

    //- Default aerosolProperties dictionary name
    static const word aerosolPropertiesName;

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        aerosolModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties
        ),
        (modelType, mesh, aerosolProperties)
    );

    // Constructors

        //- Construct from components
        aerosolModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& aerosolProperties=aerosolPropertiesName
        );

    //- Selector
    static autoPtr<aerosolModel> New
    (
        const fvMesh& mesh,
        const word& aerosolProperties=aerosolPropertiesName
    );


    //- Destructor
    virtual ~aerosolModel();


    // Member Functions

        //- Return access to the thermo package
        inline aerosolThermo& thermo();

        //- Return access to the thermo package
        inline const aerosolThermo& thermo() const;

        //- Return access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return access to the drift flux model
        inline driftFluxModel& drift();

        //- Return access to the drift flux model
        inline const driftFluxModel& drift() const;

        //- Set turbulence
        inline void setTurbulence(compressibleTurbulenceModel& turbModel);

        //- Return access to turbulence
        inline const compressibleTurbulenceModel& turbulence() const;

        //- Return access to rho
        inline const volScalarField& rho() const;

        //- Return access to phi
        inline const surfaceScalarField& phi() const;

        //- Return access to U
        inline const volVectorField& U() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct the aerosol model
        virtual void correct();

        //- Correct the underlying aerosol model
        virtual void correctModel() = 0;

        //- Solution step before the mass fraction solution
        virtual void solvePre() = 0;

        //- Solution step after the mass fraction solution
        virtual void solvePost() = 0;

        //- Right-hand side source term
        virtual tmp<fvScalarMatrix> R(const volScalarField& Y) const = 0;

        //- Heat release rate
        virtual tmp<volScalarField> Qdot() const = 0;

        //- Update properties from given dictionary
        virtual bool read();

        //- Return access to the output properties dictionary
        inline const IOdictionary& outputProperties() const;

        //- Return output properties dictionary
        inline IOdictionary& outputProperties();

        //- Return reference to multivariate convectionScheme for phi
        inline const fv::convectionScheme<scalar>& mvPhi() const;

        //- Return reference to multivariate convectionScheme for phiDrift
        inline const fv::convectionScheme<scalar>& mvPhiDrift() const;

        //- Return reference to multivariate convectionScheme for phiInertial
        inline const fv::convectionScheme<scalar>& mvPhiInertial() const;

        //- Return reference to multivariate convectionScheme for phiBrownian
        inline const fv::convectionScheme<scalar>& mvPhiBrownian() const;

        //- Return the mean diameter related to the moment p and q ratio
        virtual tmp<volScalarField> meanDiameter
        (
            const scalar p,
            const scalar q
        ) const = 0;

        //- Return the median diameter related to the moment p
        virtual tmp<volScalarField> medianDiameter(const scalar p) const = 0;

        inline scalar dMin() const
        {
            return dMin_;
        }

        inline scalar dMax() const
        {
            return dMax_;
        }

        //- Return access to the corrective mass fraction drift flux
        inline const surfaceScalarField& phiDrift() const;

        //- Return access to the mixture drift stress tensor
        inline const volSymmTensorField& tauDrift() const;

        //- Return access to the inertial dispersed mass fraction drift flux
        inline const surfaceScalarField& phiInertial() const;

        //- Return access to the Brownian dispersed mass fraction drift flux
        inline const surfaceScalarField& phiBrownian() const;

        //- Return access to the dispersed mass fraction diffusivity
        inline const volScalarField& DDisp() const;

        //- Return access to the continuous mass fraction diffusivities
        inline const PtrList<volScalarField>& DCont() const;

        //- Return access to the effective droplet number flux fields list
        inline const PtrList<surfaceScalarField>& phiEff() const;

        //- get the timescale in the cells.
        tmp<scalarField> getRDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aerosolModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
