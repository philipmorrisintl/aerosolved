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

\file condensationEvaporationModel.H
\brief Condensation/evaporation library

This library provides a condensation or evaporation rate (depending on the sign)
based on a model of choice. Access is given through the getCondRate functions.
These functions are bound to function pointers inside the aerosolModel object.
In this way the aerosolModel object can directly get the condensation rates. The
condensation rates are provided in terms of kg/s, i.e., a mass growth rate,
corresponding to \f$I(s,\mathbf{x},t)\f$ in Eq. (5.1) in \cite thesis. The
getCondRate functions are pure virtual and must be 'overloaded' by the models.

*/

#ifndef condensationEvaporationModel_H
#define condensationEvaporationModel_H

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
                     class condensationEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationEvaporationModel
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

        //- Condensation rate field
        PtrList<PtrList<volScalarField> > I_;

        //- Xi-space and Gamma per species
        PtrList<PtrList<volScalarField> > etaGamma_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        condensationEvaporationModel(const condensationEvaporationModel&);

        //- Disallow default bitwise assignment
        void operator=(const condensationEvaporationModel&);


public:

    //- Runtime type information
    TypeName("condensationEvaporationModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        condensationEvaporationModel,
        dictionary,
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        ),
        (mesh, aerosol)
    );


    // Constructors

        //- Construct from components and aerosol model
        condensationEvaporationModel(const fvMesh& mesh, aerosolModel& aerosol);


    //- Selector
    static autoPtr<condensationEvaporationModel> New
    (
        const fvMesh& mesh,
        aerosolModel& aerosol
    );


    //- Destructor
    virtual ~condensationEvaporationModel();


    // Member Functions

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return access to the thermo object
            inline aerosolModel& aerosol();

            //- Return access to the thermo object
            inline fluidThermo& thermo();


        // Evolution

            //- Return the condensation rate for a given cell, for a list of
            //- droplet sizes and a for each species (pure virtual)
            virtual List<List<scalar> > getCondRateListCell
            (
                const List<scalar>& z,
                const label jCell
            ) = 0;

            //- Return the condensation rate field for a list of droplet
            //- sizes and for each species
            PtrList<PtrList<volScalarField> >& getCondRateList
            (
                const List<scalar>& z
            );

            //- Return the condensation rate for a droplet size field and
            //- for each species
            PtrList<PtrList<volScalarField> >& getCondRateField
            (
                volScalarField& z
            );

            //- Return the condensation rate for a list of droplet size
            //- fields and for each species
            PtrList<PtrList<volScalarField> >& getCondRateFields
            (
                const PtrList<volScalarField>& z
            );

            //- Return the xi-space droplet sizes, the d-independent
            //- characteristic condensational growth velocities (per species)
            //- for the condensational growth sub problem, for a list of droplet
            //- sizes and a specific cell (pure virtual)
            virtual List<List<scalar> > getEtaGammaListCell
            (
                const List<scalar>& z,
                const label jCell
            ) = 0;

            //- Return droplet size for a given list of zetas and cell index
            virtual List<scalar> psiInv
            (
                const List<scalar>& zeta,
                const label jCell
            ) = 0;

            //- Return the xi-space droplet sizes, the d-independent
            //- characteristic condensational growth velocities (per species)
            //- for the condensational growth sub problem, for a list of droplet
            //- sizes (pure virtual)
            PtrList<PtrList<volScalarField> >& getEtaGammaList
            (
                const List<scalar>& z
            );



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

#include "condensationEvaporationModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
