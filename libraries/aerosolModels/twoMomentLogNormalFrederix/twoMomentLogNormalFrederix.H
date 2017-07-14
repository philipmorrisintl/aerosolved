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

\file twoMomentLogNormalFrederix.H
\brief A two-moment log-normal aerosol model

This is a two-moment model in which only one droplet number concentration is
solved as the zeroth moment. The other moment follows from the mass fraction
equations (in terms of a droplet size distribution based on diameter $d$, this
becomes the third moment of the size distribution). There is only one `section'
which is not a real section but contains information about the whole
distribution, e.g., the number concentration \f$M\f$, droplet velocity
\f$\mathbf{v}\f$ and number source term \f$J\f$. For each species a mass
transfer source term \f$S\f$ is computed.

The local member function `updateDistData()' updates distribution information,
such as the count median diameter (CMD) which is based on \f$M\f$, the liquid
mass fractions and the critical lower droplet diameter, and a field `Q' which
scales a unit distribution into the 'real' distribution based on the non-zero
lower critical droplet diameter. The 'update()' member function does most of the
work of this model. For each cell, a species-specific condensation rate is
computed as well as a total number concentration nucleation rate. These are then
directly used inside the solver as explicit source terms. The rate of change of
the total droplet number concentration is amended with a coagulation rate.

*/

#ifndef twoMomentLogNormalFrederix_H
#define twoMomentLogNormalFrederix_H

#include "aerosolModel.H"
#include "DataEntry.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                   Class twoMomentLogNormalFrederix Declaration
\*---------------------------------------------------------------------------*/

class twoMomentLogNormalFrederix
:
    public aerosolModel
{
    // Private data

        const dimensionSet dropletSizeDimension_;

        //- Geometric standard deviation of the particle size distribution
        scalar W_;

        //- Number of intervals to approximate the condensation integral with
        label Mint_;

        //- Critical droplet diameter
        scalar Dcrit_;

        //- List of CMD and Q field
        PtrList<volScalarField> tDistData_;

        //- max iterations to find CMD
        label maxIter_;

        //- CMD iteration tolerance
        scalar TOL_;


    // Private Member Functions

        //- Disallow copy construct
        twoMomentLogNormalFrederix(const twoMomentLogNormalFrederix&);

        //- Disallow default bitwise assignment
        void operator=(const twoMomentLogNormalFrederix&);


protected:

    // Protected data

        const fvMesh& mesh_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("twoMomentLogNormalFrederix");


    // Constructors

        //- Construct from components
        twoMomentLogNormalFrederix(const fvMesh& mesh);


    //- Destructor
    virtual ~twoMomentLogNormalFrederix();


    // Member Functions

        // Access

            //- This is a moment method
            label modType() const {return MOMENTAEROSOLMODEL;};

            //- Return the size distribution type
            label sizeDistType() const {return NOSIZEDIST;};

            //- Droplet size dimension
            const dimensionSet& dropletSizeDimension() const
            {
                return dropletSizeDimension_;
            };

            inline const scalar& W() const;

        // Functions

            //- Return the field of count mean droplet diameters
            void updateDistData();

            void updateDropletFluxes();

            void update();

            void fractionalStepInternal();
            void fractionalStepExternal();

            void checkConsistency();

            void correctSizeDistribution();

            //- Return the count mean droplet diameter
            tmp<volScalarField> dcm();

            //- Return the mass mean droplet diameter
            tmp<volScalarField> dmm();

        // I-O

            //- Read info from dictionaries
            bool read();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aerosolModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoMomentLogNormalFrederixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
