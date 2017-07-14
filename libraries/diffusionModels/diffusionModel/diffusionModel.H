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

\file diffusionModel.H
\brief Provides the per-species effective vapor diffusivity

Provides the per-species effective vapor diffusivity constant
\f$\mathbb{D}_j^v\f$ and the corrective diffusion flux based on the corrective
velocity \f$\mathbf{u}_c\f$, e.g., see Eq (2.19) in \cite thesis.

*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvmSup.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
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

        //- Diffusivity coefficients for each vapor species
        PtrList<volScalarField> DY_;

        //- Correction flux velocity
        surfaceScalarField phic_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        diffusionModel(const diffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionModel&);


public:

    //- Runtime type information
    TypeName("diffusionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        diffusionModel,
        dictionary,
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        ),
        (mesh, thermo)
    );


    // Constructors

        //- Construct from components
        diffusionModel
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Selector
    static autoPtr<diffusionModel> New
    (
        const fvMesh& mesh,
        fluidThermo& thermo
    );


    //- Destructor
    virtual ~diffusionModel();


    // Member Functions

        // Access

            //- Return binary diffusivity coefficients
            inline const PtrList<volScalarField>& DY() const;

            //- Return pointer to surface scalar field
            inline const surfaceScalarField& phic() const;

            //- Return access to the thermo object
            inline fluidThermo& thermo();

            //- Return diffusion fluxes for gases
            inline const PtrList<surfaceScalarField>& J() const;

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

#include "diffusionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

