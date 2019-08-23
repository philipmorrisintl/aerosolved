/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2019 Philip Morris International

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

\file diffusivityModel.H
\brief Base class for diffusivity model

*/

#ifndef diffusivityModel_H
#define diffusivityModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//Forward decleration
class aerosolThermo;

/*---------------------------------------------------------------------------*\
                         Class diffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class diffusivityModel
:
    public refCount
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const diffusivityModel&);


protected:

    // Protected data

        //- Name of entry
        const word name_;

        //- Reference to the thermo library
        aerosolThermo& thermo_;

        //- Species labels
        const label j_;
        const label k_;


public:

    //- Runtime type information
    TypeName("diffusivityModel")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        diffusivityModel,
        dictionary,
        (
            const word& entryName,
            const dictionary& dict,
            aerosolThermo& thermo,
            const label j,
            const label k
        ),
        (entryName, dict, thermo, j, k)
    );


    // Constructor

        //- Construct from entry name
        diffusivityModel
        (
            const word& entryName,
            aerosolThermo& thermo,
            const label j,
            const label k
        );

        //- Copy constructor
        diffusivityModel(const diffusivityModel& de);

        //- Construct and return a clone
        virtual tmp<diffusivityModel> clone() const = 0;


    //- Selector
    static autoPtr<diffusivityModel> New
    (
        const word& entryName,
        const dictionary& dict,
        aerosolThermo& thermo,
        const label j,
        const label k
    );


    //- Destructor
    virtual ~diffusivityModel();


    // Member Functions

        // Access

            //- Return the name of the entry
            const word& name() const;


        // Evaluation

            //- Return the binary diffusivity field
            virtual tmp<scalarField> D() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
