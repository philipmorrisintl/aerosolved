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

\file ManninenDrift.H
\brief Manninen's local equilibrium droplet velocity model

In the local equilibrium approximation, the droplet accelleration is taken equal
to the local fluid accelleration. This then gives an algebraic expression for
the droplet velocity, i.e.,

\f[
\mathbf{v}(s) = \mathbf{u} + \tau\left[(1-\gamma)\mathbf{g}
- \mathrm{d}_t\mathbf{u}\right],
\f]

see Eq. (6.22) in \cite thesis. This model evaluates, given a droplet size,
this algebraic expression. See Sec. 6.2.2 for more information.

The Cunningham correction factor can be included using the Cunningham parameter
(Switch) inside the driftVelocityModelParameters dict of the
driftVelocityProperties file. Gravity can be set using the g keyword.

*/

#ifndef ManninenDrift_H
#define ManninenDrift_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                       Class ManninenDrift Declaration
\*---------------------------------------------------------------------------*/

class ManninenDrift
:
    public driftVelocityModel
{
    // Private data

        //- Gravity vector
        dimensionedVector g_;

        //- Switch to enable the Cunningham correction factor
        Switch Cunningham_;


    // Private Member Functions

        //- Disallow copy construct
        ManninenDrift(const ManninenDrift&);

        //- Disallow default bitwise assignment
        void operator=(const ManninenDrift&);


protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("ManninenDrift");


    // Constructors

        //- Construct from components
        ManninenDrift(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~ManninenDrift();


    // Member Functions

        // Evolution

            void updateDropDriftVelFields();

            void updateDropDriftVelField(const volScalarField& d);

            //- Function which, given the old velocity field, diameter,
            //- drag coefficient and gravity accelleration, computes the new
            //- droplet velocity
            void updateDropletVelocity
            (
                volVectorField& V,
                const volScalarField& d,
                const volScalarField& D,
                const volVectorField& G
            );


        //  Access


        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
