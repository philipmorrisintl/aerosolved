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

\file BrownianModel.H
\brief Base class of the Brownian drift flux model

*/

#ifndef BrownianModel_H
#define BrownianModel_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "PtrList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BrownianModel Declaration
\*---------------------------------------------------------------------------*/

class BrownianModel
:
    public aerosolSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BrownianModel(const BrownianModel&);

        //- Disallow default bitwise assignment
        void operator=(const BrownianModel&);


public:

    //- Runtime type information
    TypeName("BrownianModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             BrownianModel,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict
             ),
             (aerosol, dict)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        BrownianModel
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected Brownian model
        static autoPtr<BrownianModel> New
        (
            aerosolModel& aerosol,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BrownianModel();


    // Member Functions

        //- Compute the diffusivity given a droplet size field
        virtual tmp<volScalarField> D(const volScalarField& d) const;

        //- Compute the diffusivity given a droplet size and cell
        virtual scalar D(const scalar& d, const label& celli) const = 0;

        //- Compute the diffusivity given a droplet size field and patch
        virtual scalarField D
        (
            const scalarField& d,
            const label& patchi
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
