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

\file aerosolThermo.H
\brief Main aerosolThermo class

This class provides the thermodynamic state of the two-phase
continuous-dispersed mixture. It does so, by combining two thermo objects: one
for the continuous phase and one for the dispersed phase. These two thermo
objects are mixed together using phase mixing laws which can be specified in a
case for each mixture property. The aerosolThermo class is derived from the
psiThermo class, which is a standard OpenFOAM library. This means that the
aerosolThermo class works together with other standard OpenFOAM tools, such as
turbulence models.

*/

#ifndef aerosolThermo_H
#define aerosolThermo_H

#include "psiThermo.H"
#include "rhoAerosolPhaseThermo.H"
#include "pureMixture.H"
#include "mixtureDiffusivityModel.H"
#include "phaseMixing.H"
#include "speciesTable.H"
#include "multivariateScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aerosolThermo Declaration
\*---------------------------------------------------------------------------*/

class aerosolThermo
:
    public psiThermo
{
    // Private data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Thermo-package of the continuous phase
        autoPtr<rhoAerosolPhaseThermo> thermoCont_;

        //- Thermo-package of the dispersed phase
        autoPtr<rhoAerosolPhaseThermo> thermoDisp_;

        //- Continuous mixture diffusivity models
        autoPtr<mixtureDiffusivityModel> diffusivity_;

        //- Phase mixing object
        autoPtr<phaseMixing> phaseMix_;

        //- Species table
        speciesTable species_;

        //- Active species table
        speciesTable activeSpecies_;

        //- Inactive species table
        speciesTable inactiveSpecies_;

        //- Active species map
        labelList activeSpeciesMap_;

        //- Inactive species map
        labelList inactiveSpeciesMap_;

        //- Continuous species map
        labelList contSpeciesMap_;

        //- Dispersed species map
        labelList dispSpeciesMap_;

        //- Inert specie name
        word inertSpecie_;

        //- Table of mass fraction fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

        //- Table of mass fraction fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsY_;

        //- Table of mass fraction fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fieldsZ_;


public:

    //- Runtime type information
    TypeName("aerosolThermo");


    // Constructors

        //- Construct from components
        aerosolThermo
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~aerosolThermo();


    // Member Functions

        virtual word thermoName() const
        {
            return typeName;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const inline word phaseNameCont() const
        {
            return "continuous";
        }

        const inline word phaseNameDisp() const
        {
            return "dispersed";
        }

        const rhoAerosolPhaseThermo& thermoCont() const
        {
            return thermoCont_();
        }

        const rhoAerosolPhaseThermo& thermoDisp() const
        {
            return thermoDisp_();
        }

        rhoAerosolPhaseThermo& thermoCont()
        {
            return thermoCont_();
        }

        rhoAerosolPhaseThermo& thermoDisp()
        {
            return thermoDisp_();
        }

        mixtureDiffusivityModel& diffusivity()
        {
            return diffusivity_();
        }

        phaseMixing& phaseMix()
        {
            return phaseMix_();
        }

        const speciesTable& species() const
        {
            return species_;
        }

        const speciesTable& activeSpecies() const
        {
            return activeSpecies_;
        }

        const speciesTable& inactiveSpecies() const
        {
            return inactiveSpecies_;
        }

        const labelList& activeSpeciesMap() const
        {
            return activeSpeciesMap_;
        }

        const labelList& inactiveSpeciesMap() const
        {
            return inactiveSpeciesMap_;
        }

        const speciesTable& contSpecies() const
        {
            return thermoCont().species();
        }

        const speciesTable& dispSpecies() const
        {
            return thermoDisp().species();
        }

        const labelList& contSpeciesMap() const
        {
            return contSpeciesMap_;
        }

        const labelList& dispSpeciesMap() const
        {
            return dispSpeciesMap_;
        }

        const word& inertSpecie() const
        {
            return inertSpecie_;
        }

        inline const multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        fields() const
        {
            return fields_;
        }

        inline const multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        fieldsY() const
        {
            return fieldsY_;
        }

        inline const multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        fieldsZ() const
        {
            return fieldsZ_;
        }

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Update mixture properties
        virtual void correct();

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                NotImplemented;
                return thermoCont_->he();
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                NotImplemented;
                return thermoCont_->he();
            }

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        // Fields derived from thermodynamic state variables

            //- Density [kg/m3]
            virtual tmp<volScalarField> rho() const;

            //- Patch-based density [kg/m3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Extend psiThermo with a cell-based density [kg/m3]
            virtual tmp<scalarField> rho(const labelList& cells) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Mean free path [m]
            virtual tmp<volScalarField> lambda() const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            tmp<scalarField> alphahe(const label patchi) const;


        // Mixture composition variables

            //- Const access to continuous mass fractions
            virtual inline const PtrList<volScalarField>& Y() const
            {
                return thermoCont_->composition().Y();
            }

            //- Const access to dispersed mass fractions
            virtual inline const PtrList<volScalarField>& Z() const
            {
                return thermoDisp_->composition().Y();
            }

            //- Non const access to continuous mass fractions
            virtual inline PtrList<volScalarField>& Y()
            {
                return thermoCont_->composition().Y();
            }

            //- Non const access to dispersed mass fractions
            virtual inline PtrList<volScalarField>& Z()
            {
                return thermoDisp_->composition().Y();
            }

            //- Sum of continuous mass fractions
            virtual tmp<volScalarField> sumY() const;

            //- Sum of dispersed mass fractions
            virtual tmp<volScalarField> sumZ() const;

            //- Sum of continuous mass fractions for a patch
            virtual tmp<scalarField> sumY(const label patchi) const;

            //- Sum of dispersed mass fractions for a patch
            virtual tmp<scalarField> sumZ(const label patchi) const;

            //- Sum of continuous mass fractions for a list of cells
            virtual tmp<scalarField> sumY(const labelList& cells) const;

            //- Sum of dispersed mass fractions for a list of cells
            virtual tmp<scalarField> sumZ(const labelList& cells) const;

            //- Mixture molecular weight
            virtual tmp<volScalarField> W() const;

            //- Mixture molecular weight
            virtual tmp<volScalarField> WMix() const;

            //- Mixture molecular weight for a patch
            virtual tmp<scalarField> WMix(const label patchi) const;

            //- Mixture molecular weight for a list of cells
            virtual tmp<scalarField> WMix(const labelList& cells) const;


        // Lists of species properties

            //- Saturation presures for a list of species
            virtual PtrList<scalarField> pSat
            (
                const speciesTable& species
            );

            //- Dispersed densities for a list of species
            virtual PtrList<scalarField> rhoDisp
            (
                const speciesTable& species
            ) const;

            //- Continuous densities for a list of species
            virtual PtrList<scalarField> rhoCont
            (
                const speciesTable& species
            ) const;

            //- Surface tensions for a list of species
            virtual PtrList<scalarField> sigma
            (
                const speciesTable& species
            );

    // IO

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
