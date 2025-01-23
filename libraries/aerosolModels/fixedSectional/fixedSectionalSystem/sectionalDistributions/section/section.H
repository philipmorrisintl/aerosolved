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

\file section.H
\brief Section class, containing information on a specific section

*/

#ifndef section_H
#define section_H

#include "scalar.H"
#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class section Declaration
\*---------------------------------------------------------------------------*/

class section
{
protected:

    // Protected data

        //- Section number
        label sectionNum_;

        //- Section name
        word sectionName_;

        //- Representative size
        scalar x_;

        //- Lower edge size
        scalar yl_;

        //- Upper edge size
        scalar yu_;

        //- Dimensions of the size scalar
        dimensionSet sizeDims_;

        //- Sectional number concentration field [#/kg]
        autoPtr<volScalarField> M_;

        //- Sectional inertial drift flux field [kg/s]
        surfaceScalarField phiInertial_;

        //- Sectional Brownian drift flux field [kg/s]
        surfaceScalarField phiBrownian_;

        //- Sectional diffusivity field [m^2/s]
        volScalarField D_;

        //- is the field M_ valid (either read or initialized)
        bool validM_;

public:

    // Constructors

        //- Construct from data
        section
        (
            const label,
            const word,
            const scalar,
            const scalar,
            const scalar,
            const dimensionSet,
            const fvMesh&
        );

    //- Destructor
    virtual ~section();


    // Member Functions

        // Access

            //- Section number
            inline scalar sectionNum() const
            {
                return sectionNum_;
            }

            //- Section number
            inline word sectionName() const
            {
                return sectionName_;
            }

            //- Representative size

            inline scalar& x()
            {
                return x_;
            }

            inline const scalar& x() const
            {
                return x_;
            }

            inline dimensionedScalar xd() const
            {
                return dimensionedScalar("x", sizeDims_, x_);
            }

            //- Lower representative size

            inline dimensionedScalar yld() const
            {
                return dimensionedScalar("yl", sizeDims_, yl_);
            }

            inline scalar& yl()
            {
                return yl_;
            }

            inline const scalar& yl() const
            {
                return yl_;
            }

            //- Upper representative size

            inline dimensionedScalar yud() const
            {
                return dimensionedScalar("yu", sizeDims_, yu_);
            }

            inline scalar& yu()
            {
                return yu_;
            }

            inline const scalar& yu() const
            {
                return yu_;
            }

            //- Const access to the dimensions
            inline const dimensionSet& sizeDimensions() const
            {
                return sizeDims_;
            }

            //- Return the sectional mass concentration field
            inline volScalarField& M()
            {
                return M_();
            }

            inline const volScalarField& M() const
            {
                return M_();
            }

            //- Return the sectional inertial drift flux field
            inline surfaceScalarField& phiInertial()
            {
                return phiInertial_;
            }

            inline const surfaceScalarField& phiInertial() const
            {
                return phiInertial_;
            }

            //- Return the sectional inertial drift flux field
            inline surfaceScalarField& phiBrownian()
            {
                return phiBrownian_;
            }

            inline const surfaceScalarField& phiBrownian() const
            {
                return phiBrownian_;
            }

            //- Return the sectional diffusivity field
            inline volScalarField& D()
            {
                return D_;
            }

            inline const volScalarField& D() const
            {
                return D_;
            }

            //- Compute the droplet diameter given a dispersed phase density
            inline scalar d(const scalar rho) const
            {
                return
                    Foam::pow(x_/rho*6.0/constant::mathematical::pi, 1.0/3.0);
            }

            //- Compute the droplet diameter given a dispersed phase mass
            //  density field
            inline scalarField d(const scalarField& rho) const
            {
                scalarField result(x_/rho*6.0/constant::mathematical::pi);

                for (scalar& val : result)
                {
                    val = Foam::pow(val, 1.0/3.0);
                }
                return result;
            }

            //- is M_ valid?
            inline bool validM() const
            {
                return validM_;
            }

            //- set M from the value of a patch
            void initM(const word& patch);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
