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

\file FrederixDrift.H
\brief Solution to the full unsteady Lagrangian droplet equation of motion

This model computes the evolution of the droplet velocity field \f$\mathbf{v}\f$
based on the full Lagrangian equation of motion for droplets, i.e.,

\f[
\partial_t\mathbf{v}(s) + [\mathbf{v}(s)\cdot\nabla]\mathbf{v}(s) =
-\frac{1}{\tau} [\mathbf{v}(s)-\mathbf{u}] + (1-\gamma)\mathbf{g}.
\f]

Since this is an evolution equation in time, a time discretisation must be used
requiring the need to store all \f$\mathbf{v}\f$ fields to be used in the next
time step. See Sec. 6.2.2 of \cite thesis for more information.

The SchillerNaumann keyword (Switch) inside the driftVelocityModelParameters
dict of the driftVelocityProperties controls if the Schiller-Naumann drag force
(see \cite Clift:1978aa) is used. The Cunningham keyword (Switch) adds the
Cunningham slip factor into the drag term, relevant for small droplets. Also,
maxIter controls the maximum number of iterations (i.e., \f$\mathcal{K}\f$ in
Sec. 6.3.2 of \cite thesis) and TOL sets the relative tolerance below which the
solutions for \f$\mathbf{v}\f$ must be dropped for convergence. Finally, gravity
can be set by the vector \f$\mathbf{g}\f$.

*/

#ifndef FrederixDrift_H
#define FrederixDrift_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                       Class FrederixDrift Declaration
\*---------------------------------------------------------------------------*/

class FrederixDrift
:
    public driftVelocityModel
{
    // Private data

        //- Gravity vector
        dimensionedVector g_;

        //- Solver tolerance
        scalar TOL_;

        //- Max number of solver iterations
        label maxIter_;

        //- Switch to enable Schiller-Naumann drag coefficient
        Switch SchillerNaumann_;

        //- Switch to enable the Cunningham correction factor
        Switch Cunningham_;

        //- Number of time steps before next update
        label maxSkip_;

        //- Per-section number of time steps since last update
        labelList skip_;

        //- Per-section last first initial residuals
        scalarList r0_;


    // Private Member Functions

        //- Disallow copy construct
        FrederixDrift(const FrederixDrift&);

        //- Disallow default bitwise assignment
        void operator=(const FrederixDrift&);


protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("FrederixDrift");


    // Constructors

        //- Construct from components
        FrederixDrift(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~FrederixDrift();


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
                const volVectorField& G,
                const label i
            );

            //- Particle Reynold's number
            tmp<volScalarField> Re
            (
                const volScalarField& d,
                const volVectorField& U,
                const volVectorField& V
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
