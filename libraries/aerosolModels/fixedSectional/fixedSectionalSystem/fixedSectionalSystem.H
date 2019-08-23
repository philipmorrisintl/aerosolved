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

\file fixedSectionalSystem.H
\brief Base class of the fixed sectional system

This class defines the sectional system in terms of the distribution of
sections, the sections themselves and the interpolation scheme used to
interpolate data in between representative sizes of sections. Also, it contains
additional helper functions for the sectional distribution.

*/

#ifndef fixedSectionalSystem_H
#define fixedSectionalSystem_H

#include "aerosolSubModelBase.H"
#include "sectionalDistribution.H"
#include "sectionalInterpolation.H"
#include "volFields.H"
#include "regIOobject.H"
#include "multivariateScheme.H"
#include "coalescencePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedSectionalSystem Declaration
\*---------------------------------------------------------------------------*/

class fixedSectionalSystem
:
    public aerosolSubModelBase,
    public regIOobject
{
private:

    // Private Member Data

        //- Sectional distribution
        autoPtr<sectionalDistribution> distribution_;

        //- Sectional interpolation
        autoPtr<sectionalInterpolation> interpolation_;

        //- Base sectional mass fraction field
        volScalarField M_;

        //- Table of sectional mass fraction fields for multivariate convection
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

        //- Coalescence pairs
        PtrList<coalescencePair> coalescencePairs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fixedSectionalSystem(const fixedSectionalSystem&);

        //- Disallow default bitwise assignment
        void operator=(const fixedSectionalSystem&);


public:

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        fixedSectionalSystem
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fixedSectionalSystem();


    // Access

        //- Access to the sectional distribution

        virtual inline sectionalDistribution& distribution()
        {
            return distribution_();
        }

        virtual inline const sectionalDistribution& distribution() const
        {
            return distribution_();
        }

        //- Access to the sectional interpolation scheme

        virtual inline sectionalInterpolation& interpolation()
        {
            return interpolation_();
        }

        virtual inline const sectionalInterpolation& interpolation() const
        {
            return interpolation_();
        }

        //- Access to the table of sectional mass fraction fields

        inline const multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        fields() const
        {
            return fields_;
        }

        //- Access to the coalescence pairs
        inline const PtrList<coalescencePair>& coalescencePairs() const
        {
            return coalescencePairs_;
        }


    // Member Functions

        //- Compute the sectional diameter
        virtual tmp<volScalarField> d(const label i) const;

        //- Compute the mean diameter related to the moment p and q ratio
        virtual tmp<volScalarField> meanDiameter
        (
            const scalar p,
            const scalar q
        ) const;

        //- Comptue the median diameter related to the moment p
        virtual tmp<volScalarField> medianDiameter
        (
            const scalar p
        ) const;

        //- Compute the total dispersed phase mass fraction
        virtual tmp<volScalarField> alpha() const;

        //- Compute the total number concentration
        virtual tmp<volScalarField> sumM() const;

        //- Rescale the sectional system
        virtual void rescale();

        //- Generate coalescence pairs
        virtual void generateCoalescencePairs();

        //- Write
        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
