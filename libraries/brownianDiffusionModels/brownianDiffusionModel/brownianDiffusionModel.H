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

\file brownianDiffusionModel.H
\brief Provides a diffusivity associated with each droplet section.

*/

#ifndef brownianDiffusionModel_H
#define brownianDiffusionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fluidThermo.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class brownianDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class brownianDiffusionModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to aerosol object
        aerosolModel& aerosol_;

        //- Reference to thermo object
        fluidThermo& thermo_;

        //- Coefficients
        dictionary coeffs_;

        //- Diffusivity coefficients for each droplet size
        PtrList<volScalarField> DM_;


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        brownianDiffusionModel(const brownianDiffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const brownianDiffusionModel&);


public:

    //- Runtime type information
    TypeName("brownianDiffusionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        brownianDiffusionModel,
        dictionary,
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        ),
        (mesh, aerosol)
    );


    // Constructors

        //- Construct from components
        brownianDiffusionModel
        (
            const fvMesh& mesh,
            aerosolModel& aerosol
        );


    //- Selector
    static autoPtr<brownianDiffusionModel> New
    (
        const fvMesh& mesh,
        aerosolModel& aerosol
    );


    //- Destructor
    virtual ~brownianDiffusionModel();


    // Member Functions

        // Access

            //- Return binary diffusivity coefficients
            inline const PtrList<volScalarField>& DM() const;

            //- Return access to the aerosol object
            inline aerosolModel& aerosol();

        // Thermodynamics

            //- Update effective diffusivity (pure virtual)
            virtual void update() = 0;

        // I-O

            //- Read info from dictionaries (pure virtual)
            virtual bool read() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "brownianDiffusionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

