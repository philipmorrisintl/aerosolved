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

\file driftVelocityModel.H
\brief Model which provides the droplet velocity

This library provides the droplet velocity from which the drift velocity can be
computed using:

\f[
\mathbf{v}(s) = \mathbf{u} + \mathbf{u}^\ell (s),
\f]

See Eq. (6.2) in \cite thesis. Perhaps a little confusing, the \f$\mathbf{v}\f$
fields are stored inside the aerosolModel object, as they are a part of the
physical description of the droplet phase. The chosen drift velocity model sets
these fields according to some equation of motion for the droplets.

*/

#ifndef driftVelocityModel_H
#define driftVelocityModel_H

#include "fluidThermo.H"
#include "aerosolModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      class driftVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class driftVelocityModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to aerosol model object
        aerosolModel& aerosol_;

        //- Reference to thermo object
        fluidThermo& thermo_;

        //- Coefficients
        dictionary params_;


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        driftVelocityModel(const driftVelocityModel&);

        //- Disallow default bitwise assignment
        void operator=(const driftVelocityModel&);


public:

    //- Runtime type information
    TypeName("driftVelocityModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        driftVelocityModel,
        dictionary,
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        ),
        (mesh, aerosol)
    );


    // Constructors

        //- Construct from components and aerosol model
        driftVelocityModel(const fvMesh& mesh, aerosolModel& aerosol);


    //- Selector
    static autoPtr<driftVelocityModel> New
    (
        const fvMesh& mesh,
        aerosolModel& aerosol
    );


    //- Destructor
    virtual ~driftVelocityModel();


    // Member Functions

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return access to the thermo object
            inline aerosolModel& aerosol();

            //- Return access to the thermo object
            inline fluidThermo& thermo();


        // Evolution

            //- Return the droplet velocity fields (assumes sectional)
            virtual void updateDropDriftVelFields() = 0;

            //- Return the droplet velocity field for a given field of droplet
            //- diameters
            virtual void updateDropDriftVelField(const volScalarField& d) = 0;


    // I-O

        //- Update properties from given dictionary
        virtual bool read();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "driftVelocityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
