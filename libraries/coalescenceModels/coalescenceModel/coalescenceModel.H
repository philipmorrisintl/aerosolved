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

\file coalescenceModel.H
\brief Returns a coalescence rate

Given two droplet sizes, this library returns the coalescence (a.k.a.
coagulation) rate based on some model.

*/

#ifndef coalescenceModel_H
#define coalescenceModel_H

#include "fluidThermo.H"
#include "aerosolModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvmSup.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class coalescenceModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceModel
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
        dictionary coeffs_;


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        coalescenceModel(const coalescenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceModel&);


public:

    //- Runtime type information
    TypeName("coalescenceModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        coalescenceModel,
        dictionary,
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        ),
        (mesh, aerosol)
    );


    // Constructors

        //- Construct from components and aerosol model
        coalescenceModel(const fvMesh& mesh, aerosolModel& aerosol);


    //- Selector
    static autoPtr<coalescenceModel> New
    (
        const fvMesh& mesh,
        aerosolModel& aerosol
    );


    //- Destructor
    virtual ~coalescenceModel();


    // Member Functions

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return access to the aerosol model object
            inline aerosolModel& aerosol();

            //- Return access to the thermo object
            inline fluidThermo& thermo();

        // Evolution

            //- Get a field of coalescence rates for two droplet sizes
            //- (has model-independent implementation)
            virtual tmp<volScalarField> getCoaRateField
            (
                const scalar,
                const scalar,
                const volScalarField&
            );

            //- Get the coalescence rate in a specific cell for two droplet sizes
            //- and provided Knutsen number
            virtual scalar getCoaRateCell
            (
                const scalar,
                const scalar,
                const label,
                const scalar
            ) = 0;


    // I-O

        //- Update properties from given dictionary
        virtual bool read();

        inline scalar pi() const;
        inline scalar k() const;
        inline scalar N_A() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coalescenceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
