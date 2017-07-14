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

\file fluidThermo.H
\brief The fluidThermo library base class

This library primarily describes the mixture and how the mixture behaves under
changes in temperature and pressure. The library contains the temperature field
\f$T\f$, relative pressure field \f$p_1=p-p_0\f$, density \f$\rho\f$,
compressibility ratio \f$\psi\f$, heat capacities \f$c_p\f$ and/or \f$c_v\f$ and
mixture mass fraction vectors \f$\mathbf{Y}\f$ and \f$\mathbf{Z}\f$. The library
also contains general information about all species in the 'species' dictionary.

The main member functions of the fluidThermo library are 'updateRho' and
'updatePsi'. They control the updates of the \f$\rho\f$ and \f$\psi\f$ fields,
given some model. These functions are necessary to let density 'feel' changes in
pressure, temperature and composition as expressed by \f$\mathbf{Y}\f$ and
\f$\mathbf{Z}\f$. The compressibility ratio \f$\psi\f$ is defined in the
equation equation of state as the ratio between density and pressure. An example
is given by Eq. (2.8) in \cite thesis. It can be used to write the density in
terms of a pressure, e.g., the density in the continuity equation. This then
gives an evoluation equation for pressure.

Currently there are three fluidThermo models implemented in AeroSolved:
'compressible', 'semiIncompressible' and 'incompressible'. These have a
counterpart in the aerosolEulerFoam solver, in particular in the formulation of
the pressure equation (pEqn.H). We will now discuss the three models.

*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvmSup.H"
#include "DataEntry.H"
#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    public IOdictionary
{
public:
    // Data types

        //- Vapor or liquid identifier
        enum phaseType{VAPOR, LIQUID};

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Temperature field \f$T\f$
        volScalarField T_;

        //- Pressure fluctuation field \f$p_1\f$
        volScalarField p1_;

        //- Base pressure \f$p_0\f$
        dimensionedScalar p0_;

        //- Liquid and vapor density limits
        scalar rholMin_;
        scalar rholMax_;
        scalar rhovMin_;
        scalar rhovMax_;

        //- Mixture mass density \f$\rho\f$
        volScalarField rho_;

        //- Compressibility \f$\psi\f$
        volScalarField psi_;

        //- Mixture heat capacity at constant volume \f$c_v\f$
        volScalarField CvEff_;

        //- Mixture heat capacity at constant pressure \f$c_p\f$
        volScalarField CpEff_;

        //- Vapor mole fractions \f$\mathbf{W}\f$
        PtrList<volScalarField> W_;

        //- Liquid mole fractions \f$\mathbf{X}\f$
        PtrList<volScalarField> X_;

        //- Vapor mass fractions \f$\mathbf{Y}\f$
        PtrList<volScalarField> Y_;

        //- Liquid mass fractions \f$\mathbf{Z}\f$
        PtrList<volScalarField> Z_;

        //- Universal gas constant [J/mol/K]
        dimensionedScalar R_;

        //- Molar weights
        List<scalar> M_;

        //- Critical temperatures
        List<scalar> Tc_;

        //- Phase change switch
        List<Switch> phaseChange_;

        //- Species dictionary, containing info for each species
        dictionary species_;

        //- Species dictionary, containing info for each species, for only those
        //- which may change phase
        dictionary speciesPhaseChange_;

        //- Dictionary of thermodynamical model parameters
        dictionary params_;

        //- Mass fraction tolerance for the sum of all mass fractions
        scalar massConservationTolerance_;

        //- Mass fraction correction relaxation
        scalar massConservationRelaxation_;

        //- Table of list of DataEntry objects for each species
        HashTable< PtrList<DataEntry<scalar> > > dataEntriesVapor_;
        HashTable< PtrList<DataEntry<scalar> > > dataEntriesLiquid_;

        //- Table of list of switches to check if a species property is
        //- available, for each species
        HashTable< List<Switch> > propertiesAvailableVapor_;
        HashTable< List<Switch> > propertiesAvailableLiquid_;

        //- List of diffusivityModels for each species combination
        PtrList<diffusivityModel> diffusivityModels_;

        //- Matrix to find linear diffusivity index
        SquareMatrix<label> diffIndex_;


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        fluidThermo(const fluidThermo&);

        //- Disallow default bitwise assignment
        void operator=(const fluidThermo&);

        //- Set dataEntry property for this species, if available in the
        //- species dictionary. Return true if successful.
        bool readPropertyThisSpecies
        (
            PtrList< DataEntry<scalar> >& dataEntrieList,
            const label i,
            const word propertyName,
            const phaseType phase,
            const dictionary& speciesDict
        );

        void readProperty
        (
            const word propertyName,
            const phaseType phase,
            const dictionary& speciesDict,
            const bool force
        );

        //- Return access to List of dataEntry objects to evaluate properties
        PtrList<DataEntry<scalar> >& getProperty
        (
            const word propertyName,
            const phaseType phase,
            const bool force
        );

        //- Set the list of diffusivityModels
        void getDiffusivityModels();


public:

    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        dictionary,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from components
        fluidThermo(const fvMesh& mesh);


    //- Selector
    static autoPtr<fluidThermo> New(const fvMesh& mesh);


    //- Destructor
    virtual ~fluidThermo();


    // Member Functions

        static const word dictName;

        //- Return fluidThermo object based on patch field
        static const fluidThermo& lookupThermo(const fvPatchScalarField& pf);

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return temperature
            inline volScalarField& T();

            //- Return temperature
            inline const volScalarField& T() const;

            //- Return vapor mass fraction
            inline PtrList<volScalarField>& Y();
            inline const PtrList<volScalarField>& Y() const;

            //- Return vapor mole fraction
            inline PtrList<volScalarField>& W();
            inline const PtrList<volScalarField>& W() const;

            //- Return liquid mass fraction
            inline PtrList<volScalarField>& Z();
            inline const PtrList<volScalarField>& Z() const;

            //- Return liquid mole fraction
            inline PtrList<volScalarField>& X();
            inline const PtrList<volScalarField>& X() const;

            //- Universal gas constant
            inline const dimensionedScalar& R() const;

            //- Return mixture mass density
            inline volScalarField& rho();
            inline const volScalarField& rho() const;

            //- Return mixture compressibility
            inline volScalarField& psi();

            //- Return base pressure
            inline dimensionedScalar& p0();
            inline const dimensionedScalar& p0() const;

            //- Return pressure fluctuation
            inline volScalarField& p1();
            inline const volScalarField& p1() const;

            //- Return mixture heat capacity at constant volume
            inline const volScalarField& CvEff() const;

            //- Return mixture heat capacity at constant pressure
            inline const volScalarField& CpEff() const;

            //- Return species dictionary
            inline const dictionary& species() const;

            //- Return species dictionary
            inline const dictionary& speciesPhaseChange() const;

            //- Return the number of species
            inline label nSpecies() const;

            //- Return the number of phase changing species
            inline label nSpeciesPhaseChange() const;

            //- Return list of mole weights
            inline List<scalar> M() const;

            //- Return list of critical temperatures
            inline List<scalar> Tc() const;

            //- Return list of phase change switched
            inline List<Switch> phaseChange() const;

            //- Check if the specified property for the given species in the
            //- given phase is available
            bool propertyFound(const label i, const word propertyName, const phaseType phase);

            //- Return access to List of dataEntry objects to evaluate properties,
            //- if present. Else return an empty PtrList.
            PtrList<DataEntry<scalar> >& getPropertyIfPresent
            (
                const word propertyName,
                const phaseType phase
            )
            {
                return getProperty(propertyName, phase, false);
            }

            //- Return access to List of dataEntry objects to evaluate properties,
            //- if present. Else and error is thrown.
            PtrList<DataEntry<scalar> >& getProperty
            (
                const word propertyName,
                const phaseType phase
            )
            {
                return getProperty(propertyName, phase, true);
            }

            //- Matrix to find linear diffusivity index
            inline SquareMatrix<label> diffIndex() const;



        // Thermodynamics

            //- Correct species mass fractions: no negativity and they should
            //- sum up to one
            virtual void rescale
            (
                PtrList<volScalarField>& Y,
                PtrList<volScalarField>& Z,
                PtrList<volScalarField>& M,
                bool printInfo = false,
                bool clipToZero = true
            );

            //- Based on species mass fractions, update mole fractions
            virtual void updateMoleFractions();

            //- Update density, based on psi and pressure (pure virtual)
            virtual void updateRho() = 0;

            //- Update compressibility (pure virtual)
            virtual void updatePsi() = 0;

            //- Update mixture heat capacity at constant volume (pure virtual)
            virtual void updateCvEff() = 0;

            //- Update mixture heat capacity at constant pressure (pure virtual)
            virtual void updateCpEff() = 0;

            //- Total vapor mass fraction in mixture
            virtual tmp<volScalarField> Ytot();

            //- Total liquid mass fraction in mixture
            virtual tmp<volScalarField> Ztot();

            //- Current liquid density
            virtual tmp<volScalarField> rhoLiquid() = 0;

            //- Current vapor density
            virtual tmp<volScalarField> rhoVapor() = 0;

            //- Current vapor mean molecular mass
            virtual tmp<volScalarField> mVapor();

            //- Vapor density (pure virtual)
            virtual tmp<volScalarField> rhoY() = 0;

            //- Liquid density (pure virtual)
            virtual tmp<volScalarField> rhoZ() = 0;

            //- Vapor compressibility (pure virtual)
            virtual tmp<volScalarField> psiY() = 0;

            //- Current mixture density (pure virtual)
            virtual tmp<volScalarField> rhoMix() = 0;

            //- Get binary diffusivity for species a and b, based on T and p
            scalar getDiffusivity
            (
                const label a,
                const label b,
                const scalar T,
                const scalar p
            );

            //- Get field of binary diffusivities for species a and b, based on
            //- T and p
            tmp<scalarField> getDiffusivity
            (
                const label a,
                const label b,
                const scalarField& T,
                const scalarField& p
            );

            //- Get volScalarField of binary diffusivities for species a and b,
            //- based on T and p
            tmp<volScalarField> getDiffusivity
            (
                const label a,
                const label b,
                const volScalarField& T,
                const volScalarField& p
            );

            //- Clip liquid density to specified limits
            void limitVaporDensity(surfaceScalarField& rho);
            void limitLiquidDensity(surfaceScalarField& rho);
            void limitVaporDensity(volScalarField& rho);
            void limitLiquidDensity(volScalarField& rho);
            void limitVaporDensity(scalarField& rho);
            void limitLiquidDensity(scalarField& rho);


    // I-O

        //- Read info from dictionaries
        virtual bool read();

        //- Read property for specific phase and species if present
        void readPropertyIfPresent(const word propertyName, const phaseType phase, const dictionary& speciesDict)
        {
            readProperty(propertyName, phase, speciesDict, false);
        }

        //- Read property for both phases and specific species if present
        void readPropertyBothIfPresent(const word propertyName, const dictionary& speciesDict)
        {
            readProperty(propertyName, VAPOR, speciesDict, false);
            readProperty(propertyName, LIQUID, speciesDict, false);
        }

        //- Force read of property for specific phase and species
        void readProperty(const word propertyName, const phaseType phase, const dictionary& speciesDict)
        {
            readProperty(propertyName, phase, speciesDict, true);
        }

        //- Force read of property for both phases and specific species
        void readPropertyBoth(const word propertyName, const dictionary& speciesDict)
        {
            readProperty(propertyName, VAPOR, speciesDict, true);
            readProperty(propertyName, LIQUID, speciesDict, true);
        }

        //- Read property for specific phase if present, for all species
        void readPropertyIfPresent(const word propertyName, const phaseType phase)
        {
            readProperty(propertyName, phase, species_, false);
        }

        //- Read property for both phases if present, for all species
        void readPropertyBothIfPresent(const word propertyName)
        {
            readProperty(propertyName, VAPOR, species_, false);
            readProperty(propertyName, LIQUID, species_, false);
        }

        //- Force read of property for specific phase, for all species
        void readProperty(const word propertyName, const phaseType phase)
        {
            readProperty(propertyName, phase, species_, true);
        }

        //- Force read of property for both phases, for all species
        void readPropertyBoth(const word propertyName)
        {
            readProperty(propertyName, VAPOR, species_, true);
            readProperty(propertyName, LIQUID, species_, true);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

