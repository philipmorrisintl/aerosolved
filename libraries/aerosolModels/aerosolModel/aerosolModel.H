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

\file aerosolModel.H
\brief Aerosol models

This library is concerned with all information regarding that what makes a
mixture an aerosol; it contains a list of \f$M_i\f$ number concentration fields,
a list of \f$\mathbf{v}_i\f$ velocity fields, aerosol source terms and general
information about the description of the disperse character of an aerosol. Also
a number of general member functions are included, such as the
`updateSizeDistribution' routine which sets the distribution of sections given a
certain type of distribution (e.g., `none', `linear', `logarithmic', `geometric'
or `list') or the `limitWallFlux' function which clips a flux field to
positivity at all patches which are of type `wall'.

The aerosolModels library contains a pointer of type `fluidThermo'. Upon
construction of the aerosolModels library a fluidThermo library is also
automatically created and the pointer is set to point to this fluidThermo
library. This means that the user should not create a new separate fluidThermo
object next to the aerosolModels object, but should rather use the fluidThermo
object from the aerosol model. This is accessible via the `thermo()' member
function of the aerosolModels library.

*/

#ifndef aerosolModel_H
#define aerosolModel_H

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvm.H"
#include "fluidThermo.H"
#include "fvMatrices.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Model type identifier
enum aerosolModelType
{
    MOMENTAEROSOLMODEL = 0,
    SECTIONALAEROSOLMODEL = 1
};

//- Droplet size distribution type
enum aerosolSizeDistributionType
{
    NOSIZEDIST = 0,
    LINEARSIZEDIST = 1,
    LOGARITHMICSIZEDIST = 2,
    GEOMETRICSIZEDIST = 3,
    LIST = 4
};

enum aerosolSizePositionType
{
    NOSIZEPOS = 0,
    CENTERSIZEPOS = 1,
    INTERFACESIZEPOS = 2
};

enum distMethod
{
    TWOMOMENT = 0,
    FOURMOMENT = 1,
    HYBRID = 2
};

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class aerosolModel Declaration
\*---------------------------------------------------------------------------*/

class aerosolModel
:
    public IOdictionary
{
public:

    // Data types

        static const NamedEnum<aerosolModelType, 2> modelTypeNames;
        static const NamedEnum<aerosolSizeDistributionType, 5> sizeDistributionTypeNames;
        static const NamedEnum<aerosolSizePositionType, 3> sizePositionTypeNames;
        static const NamedEnum<distMethod, 3> distMethodNames;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Pointer to fluid thermo model
        autoPtr<fluidThermo> thermoPtr_;

        //- Droplet number concentrations per unit mass \f$M_i\f$
        PtrList<volScalarField> M_;

        //- Pointer to droplet velocity fields \f$\mathbf{v}_i\f$
        PtrList<volVectorField> V_;

        //- Total mass transfer rate for each species
        PtrList<volScalarField> S_;

        //- Total droplet number concentration source term
        PtrList<volScalarField> J_;

        //- List of droplet fluxes
        PtrList<surfaceScalarField> phid_;

        //- Enthalpy of vaporization source term
        volScalarField HvapS_;

        //- List of section boundary positions \f$y_i\f$
        List<scalar> y_;

        //- List of representative droplet sizes \f$x_i\f$, for each section
        List<scalar> x_;

        //- Number of sections \f$\mathcal{P}\f$
        label P_;

        //- Logarithmic size distribution power (see equation (18) in
        //- Lister:1995aa, although, in contrast with this work, here q_ is for
        //- generality a scalar, not an integer. Also, to scale the
        //- distribution between yMin and yMax, the factor 2 is replaced with
        //- a, which is computed using a = (yMax/yMin)^(q/P))
        scalar q_;

        //- Lower limit of the size distribution
        scalar yMin_;

        //- Upper limit of the size distribution
        scalar yMax_;

        //- Type of droplet size distribution
        aerosolSizeDistributionType sizeDistType_;

        //- Dictionary of aerosol model parameters
        dictionary params_;

        //- Switches
        Switch doDrift_;
        Switch doCoa_;
        Switch doNuc_;
        Switch doCond_;
        Switch doCorrSizeDist_;
        Switch doMonitors_;

        //- Monitor field hashed pointer table
        HashPtrTable<volScalarField> scalarMonitorPtrs_;
        HashPtrTable<volVectorField> vectorMonitorPtrs_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        aerosolModel(const aerosolModel&);

        //- Disallow default bitwise assignment
        void operator=(const aerosolModel&);


public:

    //- Runtime type information
    TypeName("aerosolModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        aerosolModel,
        dictionary,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from components
        aerosolModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<aerosolModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~aerosolModel();

    //- Function pointers

        //- Return the drop velocity fields for a given list of droplet sizes
        boost::function< void () > updateDropDriftVelFields;

        //- Return the drop velocity field for a given field of droplet diameters
        boost::function< void
        (
            const volScalarField&
        ) > updateDropDriftVelField;

        //- Return the drop diffusion velocity fields for a given list of
        //- droplet sizes
        boost::function< PtrList<surfaceScalarField>& () > getDropDiffFluxFields;

        //- Return the drop diffusion velocity field for a given field of
        //- droplet diameters
        boost::function< PtrList<surfaceScalarField>&
        (
            const volScalarField&
        ) > getDropDiffFluxField;

        //- Return the condensation rate for a given cell, for a list of droplet
        //- diameters and a for each species
        boost::function< List<List<scalar> >
        (
            const List<scalar>&,
            const label
        ) > getCondRateListCell;

        //- Return the condensation rate field for a list of droplet diameters
        //- and for each species
        boost::function< PtrList<PtrList<volScalarField> >&
        (
            const List<scalar>&
        ) > getCondRateList;

        //- Return the condensation rate for a list of droplet diameter fields
        //- and for each species
        boost::function< PtrList<PtrList<volScalarField> >&
        (
            const PtrList<volScalarField>&
        ) > getCondRateFields;

        //- Return the condensation rate for a droplet diameter field and for
        //- each species
        boost::function< PtrList<PtrList<volScalarField> >&
        (
            volScalarField&
        ) > getCondRateField;

        //- Return the \f$\eta\f$-space droplet sizes and the \f$d\f$-independent
        //- characteristic condensational growth velocities for each species,
        //- for the condensational growth sub problem for a list of droplet
        //- diameters and given cell
        boost::function< List<List<scalar> >
        (
            const List<scalar>&,
            const label
        ) > getEtaGammaListCell;

        //- Return the \f$\eta\f$-space droplet sizes and the \f$d\f$-independent
        //- characteristic condensational growth velocities for each species,
        //- for the condensational growth sub problem for a list of droplet
        //- diameters
        boost::function< PtrList<PtrList<volScalarField> >&
        (
            const List<scalar>&
        ) > getEtaGammaList;

        //- Return droplet diameter for a given list of zetas and cell index
        boost::function< List<scalar>
        (
            const List<scalar>&,
            const label
        ) > psiInv;

        //- Return a coalescence rate field for a given combination of droplet
        //- sizes
        boost::function< tmp<volScalarField>
        (
            const scalar,
            const scalar,
            const volScalarField&
        ) > getCoaRateField;

        //- Return the coalescence rate for a given combination of droplet sizes
        //- and cell
        boost::function< scalar
        (
            const scalar,
            const scalar,
            const label,
            const scalar
        ) > getCoaRateCell;

        //- Return the nucleation droplet and mass rates and critical nucleation
        //- cluster diameter
        boost::function< PtrList<volScalarField>& () > getNucFields;


    // Member Functions

        static const word dictName;

        // Access

            //- Return access to the thermo package
            virtual fluidThermo& thermo();

            //- Return const access to the thermo package
            virtual const fluidThermo& thermo() const;

            //- Return access to the mesh
            const fvMesh& mesh() const;

            //- Return section droplet number concentration per unit mass
            //- \f$M_i\f$
            inline PtrList<volScalarField>& M();
            inline const PtrList<volScalarField>& M() const;

            //- Return total mass transfer rate for each species
            inline PtrList<volScalarField>& S();
            inline const PtrList<volScalarField>& S() const;

            //- Return total nucleation rates
            inline PtrList<volScalarField>& J();
            inline const PtrList<volScalarField>& J() const;

            //- Return droplet flux fields
            inline PtrList<surfaceScalarField>& phid();
            inline const PtrList<surfaceScalarField>& phid() const;

            //- Return section droplet velocity \f$\mathbf{v}_i\f$
            inline PtrList<volVectorField>& V();
            inline const PtrList<volVectorField>& V() const;

            //- Return heat of vaporization
            inline volScalarField& HvapS();
            inline const volScalarField& HvapS() const;

            //- Return the type of model (pure virtual)
            virtual label modType() const = 0;

            //- Return the size distribution type (pure virtual)
            virtual label sizeDistType() const = 0;

            //- Return a reference to the dimension of the droplet size
            //- (e.g., mass, radius, volume), if applicable (pure virtual)
            virtual const dimensionSet& dropletSizeDimension() const = 0;

            //- Return a reference to the representative droplet sizes, for each
            //- section
            inline const List<scalar>& x() const;

            //- Return a reference to the section boundary positions
            inline const List<scalar>& y() const;

            //- Return the number of sections
            inline label P() const;

            //- Return the logarithmic size distribution power
            inline label q() const;

            //- Return parameter dictionary
            inline dictionary& params();
            inline const dictionary& params() const;

            //- Return const access to doMonitors_
            inline const Switch& doMonitors() const;

            //- Return const access to doDrift
            inline const Switch& doDrift() const;

            //- Return const access to monitor fields
            inline const HashPtrTable<volScalarField>& scalarMonitors() const;
            inline const HashPtrTable<volVectorField>& vectorMonitors() const;


        // Functions

            //- Update the size distribution. The default one accepts 'none',
            //- 'linear' and 'sectional' as aerosolSizeDistributionType. May be
            //- overwritten by model.
            virtual void updateSizeDistribution();

            //- Update
            virtual void update() = 0;

            //- Fractional step
            virtual void fractionalStepInternal() = 0;
            virtual void fractionalStepExternal() = 0;

            //- Check consistency of the model
            virtual void checkConsistency() = 0;

            //- Correct size distribution
            virtual void correctSizeDistribution() = 0;

            //- Compute count mean diameter (pure virtual)
            virtual tmp<volScalarField> dcm() = 0;

            //- Compute mass mean diameter (pure virtual)
            virtual tmp<volScalarField> dmm() = 0;

            //- Limit wall fluxes
            virtual void limitWallFlux(surfaceScalarField& phi) const;

            //- Set monitor fields
            void setMonitors();

            //- Clear monitor fields
            void clearMonitors();

            //- Clear scalar monitor fields
            void clearScalarMonitors();

            //- Clear vector monitor fields
            void clearVectorMonitors();

            //- Update heat of vaporization term
            void updateHvapS();

            //- Check if the updateDropDriftVelFields function pointer is set
            inline bool checkUpdateDropDriftVelFields() const;

            //- Check if the updateDropDriftVelField function pointer is set
            inline bool checkUpdateDropDriftVelField() const;

            //- Check if the getDropDiffFluxFields function pointer is set
            inline bool checkGetDropDiffFluxFields() const;

            //- Check if the getDropDiffFluxField function pointer is set
            inline bool checkGetDropDiffFluxField() const;

            //- Check if the getCondRateListCell function pointer is set
            inline bool checkGetCondRateListCell() const;

            //- Check if the getCondRateList function pointer is set
            inline bool checkGetCondRateList() const;

            //- Check if the getCondRateField function pointer is set
            inline bool checkGetCondRateField() const;

            //- Check if the getCondRateFields function pointer is set
            inline bool checkGetCondRateFields() const;

            //- Check if the getEtaGammaListCell function pointer is set
            inline bool checkGetEtaGammaListCell() const;

            //- Check if the getEtaGammaList function pointer is set
            inline bool checkGetEtaGammaList() const;

            //- Check if the getCoaRateField function pointer is set
            inline bool checkGetCoaRateField() const;

            //- Check if the getCoaRateCell function pointer is set
            inline bool checkGetCoaRateCell() const;

            //- Check if the getNucFields function pointer is set
            inline bool checkGetNucFields() const;

            //- Check if the updateDropDriftVelFields function pointer is set, else exit
            inline void checkUpdateDropDriftVelFieldsOrExit() const;

            //- Check if the updateDropDriftVelField function pointer is set, else exit
            inline void checkUpdateDropDriftVelFieldOrExit() const;

            //- Check if the getDropDiffFluxFields function pointer is set, else exit
            inline void checkGetDropDiffFluxFieldsOrExit() const;

            //- Check if the getDropDiffFluxField function pointer is set, else exit
            inline void checkGetDropDiffFluxFieldOrExit() const;

            //- Check if the getCondRateListCell function pointer is set, else exit
            inline void checkGetCondRateListCellOrExit() const;

            //- Check if the getCondRateList function pointer is set, else exit
            inline void checkGetCondRateListOrExit() const;

            //- Check if the getCondRateField function pointer is set, else exit
            inline void checkGetCondRateFieldOrExit() const;

            //- Check if the getCondRateFields function pointer is set, else exit
            inline void checkGetCondRateFieldsOrExit() const;

            //- Check if the getEtaGammaListCell function pointer is set, else exit
            inline void checkGetEtaGammaListCellOrExit() const;

            //- Check if the getEtaGammaList function pointer is set, else exit
            inline void checkGetEtaGammaListOrExit() const;

            //- Check if the getCoaRateField function pointer is set, else exit
            inline void checkGetCoaRateFieldOrExit() const;

            //- Check if the getCoaRateCell function pointer is set, else exit
            inline void checkGetCoaRateCellOrExit() const;

            //- Check if the getNucFields function pointer is set, else exit
            inline void checkGetNucFieldsOrExit() const;


    // I-O

        //- Read info from dictionaries
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aerosolModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

