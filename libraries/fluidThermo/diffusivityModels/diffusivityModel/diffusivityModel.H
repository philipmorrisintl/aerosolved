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

\file diffusivityModel.H
\brief Class which provides a model for binary diffusivity

Provides an interface to return a binary diffusivity value for a pair of two
species, independent of the mixture composition. This library is part of the
fluidThermo library, and \f$\mathcal{N}\times(\mathcal{N}-1)\f$ pairs of species
may be formed such that an equal number of diffusivityModel objects are created.

*/

#ifndef diffusivityModel_H
#define diffusivityModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "Field.H"
#include "scalar.H"
#include "refCount.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

        //- Copy of species dict
        const dictionary species_;

        //- Index in species dict to first and second species
        const label a_;
        const label b_;


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
            const dictionary species,
            const label a,
            const label b
        ),
        (entryName, dict, species, a, b)
    );


    // Constructor

        //- Construct from entry name
        diffusivityModel
        (
            const word& entryName,
            const dictionary species,
            const label a,
            const label b
        );

        //- Copy constructor
        diffusivityModel(const diffusivityModel& de);

        //- Construct and return a clone
        virtual tmp<diffusivityModel> clone() const
        {
            return tmp<diffusivityModel>(new diffusivityModel(*this));
        }


    //- Selector
    static autoPtr<diffusivityModel> New
    (
        const word& entryName,
        const dictionary& dict,
        const dictionary species,
        const label a,
        const label b
    );


    //- Destructor
    virtual ~diffusivityModel();


    // Member Functions

        // Access

            //- Return the name of the entry
            const word& name() const;


        // Evaluation

            //- Return diffusivity as a function of temperature and pressure
            virtual scalar value(const scalar T, const scalar p) const;

            //- Return diffusivity as a function of temperature and pressure
            virtual tmp<scalarField> value
            (
                const scalarField& T,
                const scalarField& p
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
