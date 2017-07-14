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

\file sectionalFrederix.H
\brief The Frederix sectional aerosol model

This is the main sectional aerosol model. The constructor of this model sets a
number of \f$\mathcal{P}\f$ (as specified in the 'aerosolProperties' dictionary
in the 'constant' case directory) \f$M_i\f$ number concentration fields. These
fields are read from disk if the 'M.0' field exists in the current time
directory, or freshly created from the 'M' field otherwise. The same strategy is
adopted for the creation of the \f$\mathbf{v}_i\f$ droplet velocity fields.

In the 'update()' member function a number of models are updated and the
required rates or coefficients are obtained from these models. These are:
nucleation rate and size, droplet drift velocities and condensation or
evaporation rates.

In the 'fractionalStepInternal()' member function Eq. (3.45) from \cite thesis
is solved. The right-hand side source term \f$\mathcal{J}_{\mathbf{X}}\f$
contains the contributions of condensation or evaporation and nucleation. The
algorithm used to solve Eq. (3.45) from \cite thesis is the CBSM, which is
documented and published in \cite Frederix:2016ae.

In the 'fractionalStepExternal()' member function a third fractional step is
solved for which the contributions have no counterpart in the mass fraction
equations and therewith have no coupling to pressure, velocity or temperature.
In this case this regards coagulation.

*/

#ifndef sectionalFrederix_H
#define sectionalFrederix_H

#include "aerosolModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                            Class sectionalFrederix Declaration
\*---------------------------------------------------------------------------*/

class sectionalFrederix
:
    public aerosolModel
{
    // Private data

        const dimensionSet dropletSizeDimension_;

        //- Pointer to eta Gamma S fields
        PtrList< PtrList<volScalarField> > *condRateListPtr_;

        //- Pointer to nucleation fields
        PtrList<volScalarField> *SJDnucPtr_;

        //- Old N fields, required for condensational growth
        PtrList<volScalarField> M0_;

        //- Old Z fields, required for condensational growth
        PtrList<volScalarField> Z0_;

        //- The method to use for distributing droplets
        distMethod distMethod_;

        //- maximum phi used in the hybrid method
        scalar phi_;

        //- Switch to check if pre-computed coalescence lists are available
        Switch preparedCoa_;

        //- Pre-computed list of coalescence zStars
        List<scalar> zStarCoa_;

        //- Pre-computed list of coalescence combinations
        List< List<label> > ijCoa_;

        //- Pre-computed list of coalescence new indices
        List<label> kCoa_;

        //- Pre-computed list of coalescence weights
        List< List<scalar> > weightsCoa_;

        //- Scalar to track the mass density of droplets which were not treated
        //- because they grew beyond the last representative size
        scalar domainDefect_;

        //- Switch to solve in zeta-space or not
        Switch zeta_;


    // Private Member Functions

        //- Disallow copy construct
        sectionalFrederix(const sectionalFrederix&);

        //- Disallow default bitwise assignment
        void operator=(const sectionalFrederix&);


protected:

    // Protected data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("sectionalFrederix");


    // Constructors

        //- Construct from components
        sectionalFrederix(const fvMesh& mesh);


    //- Destructor
    virtual ~sectionalFrederix();


    // Member Functions

        // Access

            //- This is a sectional method
            label modType() const {return SECTIONALAEROSOLMODEL;};

            //- Return the size distribution type
            label sizeDistType() const {return sizeDistType_;};

            //- Droplet size dimension
            const dimensionSet& dropletSizeDimension() const
            {
                return dropletSizeDimension_;
            };

        // Functions

            void update();

            void updateDropletFluxes();

            void fractionalStepInternal();
            void fractionalStepExternal();

            void checkConsistency();

            void correctSizeDistribution();

            tmp<volScalarField> dcm();

            tmp<volScalarField> dmm();

            //- Compute two-moment weights
            inline List<scalar> weights
            (
                scalar& x1,
                scalar& x2,
                scalar& xs
            );

            //- Compute four-moment weights
            inline List<scalar> weights
            (
                scalar& x1,
                scalar& x2,
                scalar& x3,
                scalar& x4,
                scalar& xs
            );

            //- Redistribute droplet while preserving two moments (see Sec.
            //- 4.3.1 in \cite thesis)
            inline void twoMoment
            (
                label i,
                label j,
                scalar d,
                scalar G,
                scalar phi = 1.0
            );

            //- Redistribute droplet while preserving four moments (see Sec.
            //- 4.3.1 in \cite thesis)
            inline void fourMoment
            (
                label i,
                label j,
                scalar d,
                scalar G,
                scalar phi = 1.0
            );

            //- Store the current M fields into M0_
            void storeM0();

            //- Store the current Z fields into Z0_
            void storeZ0();

            //- Get a list of indices of the section interfaces closest to the
            //- list of ds while being smaller than ds
            inline List<label> xLowerIndex(List<scalar>& ds) const;

            //- Get the index of a section interface closent to ds while being
            //- smaller than ds
            inline label xLowerIndex(scalar ds) const;

            //- Prepare coagulation (computes connectivity between sections)
            void prepareCoa();

        // I-O

            //- Read info from dictionaries
            bool read();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aerosolModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sectionalFrederixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
