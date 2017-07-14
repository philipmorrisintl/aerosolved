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

\file nucleationModel.H
\brief Provides the nucleation rate

This library provides a nucleation rate given a supersaturated state. The pure
virtual getNucFields function is registered inside the aerosolModel object using
bind such that the aerosolModel object has direct access to the model describing
the nucleation rate. The getNucFields function returns a list of volScalarFields
containing, in this order:

- A number of \f$\mathcal{N}\f$ mass transfer rates associated with the
  nucleation, for each species \f$j\f$.

- A nucleation rate indicating the number of droplets formed per unit of time

- The critical droplet size \f$z^{\mathrm{nuc}}\f$ at which nucleation takes
  place, see Eq. (5.1) in \cite thesis.

*/

#ifndef nucleationModel_H
#define nucleationModel_H

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
                     class nucleationModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to thermo object
        fluidThermo& thermo_;

        //- Reference to aerosol object
        aerosolModel& aerosol_;

        //- Coefficients
        dictionary coeffs_;

        //- Nucleation rate
        PtrList<volScalarField> SJDnuc_;


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        nucleationModel(const nucleationModel&);

        //- Disallow default bitwise assignment
        void operator=(const nucleationModel&);


public:

    //- Runtime type information
    TypeName("nucleationModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationModel,
        dictionary,
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        ),
        (mesh, aerosol)
    );


    // Constructors

        //- Construct from components
        nucleationModel(const fvMesh& mesh, aerosolModel& aerosol);


    //- Selector
    static autoPtr<nucleationModel> New
    (
        const fvMesh& mesh,
        aerosolModel& aerosol
    );


    //- Destructor
    virtual ~nucleationModel();


    // Member Functions

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return access to the thermo object
            inline fluidThermo& thermo();

            //- Return access to the aerosol object
            inline aerosolModel& aerosol();

            //- Return nucleation fields
            inline const PtrList<volScalarField>& SJDnuc() const;


        // Evolution

            //- Return the nucleation droplet and mass rates, and cirtical
            //- droplet diameter (pure virtual)
            virtual PtrList<volScalarField>& getNucFields() = 0;

        // Constants

            inline scalar k() const;
            inline scalar N_A() const;


    // I-O

        //- Update properties from given dictionary
        virtual bool read();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nucleationModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
