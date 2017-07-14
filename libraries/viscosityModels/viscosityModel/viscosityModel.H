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

\file viscosityModel.H
\brief Provides the mixture viscosity

The viscosity can be updated using the update() function and then be called by
the muEff() function.

*/

#ifndef viscosityModel_H
#define viscosityModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvmSup.H"
#include "fluidThermo.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to thermo object
        fluidThermo& thermo_;

        //- Coefficients
        dictionary coeffs_;

        //- Effective viscosity field
        volScalarField muEff_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        viscosityModel(const viscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&);


public:

    //- Runtime type information
    TypeName("viscosityModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        viscosityModel,
        dictionary,
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        ),
        (mesh, thermo)
    );


    // Constructors

        //- Construct from components
        viscosityModel
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Selector
    static autoPtr<viscosityModel> New
    (
        const fvMesh& mesh,
        fluidThermo& thermo
    );


    //- Destructor
    virtual ~viscosityModel();


    // Member Functions

        // Access

            //- Return effective viscosity
            inline const volScalarField& muEff() const;

            //- Return access to the thermo object
            inline fluidThermo& thermo();

        // Thermodynamics

            //- Update effective viscosity (pure virtual)
            virtual void update() = 0;

        // I-O

            //- Read info from dictionaries (pure virtual)
            virtual bool read() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "viscosityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
