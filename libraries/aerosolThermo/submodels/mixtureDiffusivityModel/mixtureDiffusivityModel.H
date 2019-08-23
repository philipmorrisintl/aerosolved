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

\file mixtureDiffusivityModel.H
\brief Class to compute effective mixture diffusivities

Provides the effective diffusivity of a species \f$j\f$ in the current mixture.

*/

#ifndef mixtureDiffusivityModel_H
#define mixtureDiffusivityModel_H

#include "fvMesh.H"
#include "diffusivityModel.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//Forward decleration
class aerosolThermo;

/*---------------------------------------------------------------------------*\
                       Class mixtureDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class mixtureDiffusivityModel
{
private:

    // Private Member Objects

        //- Const reference to the thermo model
        aerosolThermo& thermo_;

        //- Pointer list of diffusivity models
        PtrList<diffusivityModel> diffusivities_;

        //- Matrix to find linear diffusivity index
        SquareMatrix<label> indices_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixtureDiffusivityModel(const mixtureDiffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const mixtureDiffusivityModel&);

        //- Get the diffusivity models
        void getDiffusivityModels();


public:

    // Constructors

        //- Construct from type name, dictionary and thermo model
        mixtureDiffusivityModel
        (
            aerosolThermo& thermo
        );


    //- Destructor
    virtual ~mixtureDiffusivityModel();


    // Member Functions

        // Evolution

            //- Compute the effective continuous species diffusivity
            virtual tmp<scalarField> Deff(const label& j) const;

            //- Compute a list of effective continuous species diffusivity
            virtual PtrList<scalarField> Deff() const;

            //- Get the binary diffusivity pair by name
            virtual tmp<scalarField> D
            (
                const word& jName,
                const word& kName
            ) const;

            //- Get the binary diffusivity pair by index
            virtual tmp<scalarField> D
            (
                const label& j,
                const label& k
            ) const
            {
                return diffusivities_[indices_[j][k]].D();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
