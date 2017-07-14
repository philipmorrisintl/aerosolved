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

\file conductivityModel.H
\brief Provides the mixture heat conductivity

The effective conductivity can be updated using the update() function inside the
model. The effective conductivity can be accessed by the kEff() member function.

*/

#ifndef conductivityModel_H
#define conductivityModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      class conductivityModel Declaration
\*---------------------------------------------------------------------------*/

class conductivityModel
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

        //- Effective conductivity field
        volScalarField kEff_;

private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        conductivityModel(const conductivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityModel&);


public:

    //- Runtime type information
    TypeName("conductivityModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        conductivityModel,
        dictionary,
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        ),
        (mesh, thermo)
    );


    // Constructors

        //- Construct from components
        conductivityModel
        (
            const fvMesh& mesh,
            fluidThermo& thermo
        );


    //- Selector
    static autoPtr<conductivityModel> New
    (
        const fvMesh& mesh,
        fluidThermo& thermo
    );


    //- Destructor
    virtual ~conductivityModel();


    // Member Functions

        // Access

            //- Return effective conductivity
            inline const volScalarField& kEff() const;

            //- Return access to the thermo object
            inline fluidThermo& thermo();

        // Thermodynamics

            //- Update effective conductivity (pure virtual)
            virtual void update() = 0;

        // I-O

            //- Read info from dictionaries (pure virtual)
            virtual bool read() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conductivityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

