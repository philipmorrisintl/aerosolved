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

\file multiSpeciesFrederix.H
\brief Multi-species classical nucleation theory model

A multi-species nucleation model based on the work of \cite Winkelmann:2017aa
and further documented in \cite Frederix:2017ab. This model contains a
consistency fix with the 'saturatedVaporMixture' BC. Also, in the multi-species
formulation the mole fractions of the critical cluster now always add up to
unity (in contrast to the multiSpeciesCNT model).

*/

#ifndef multiSpeciesFrederix_H
#define multiSpeciesFrederix_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                     Class multiSpeciesFrederix Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesFrederix
:
    public nucleationModel
{
    // Private data

        //- Tolerance for detecting supersaturation
        scalar TOL_;

    // Private Member Functions

        //- Disallow copy construct
        multiSpeciesFrederix(const multiSpeciesFrederix&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesFrederix&);


protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("multiSpeciesFrederix");


    // Constructors

        //- Construct from components
        multiSpeciesFrederix(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~multiSpeciesFrederix();


    // Member Functions

        // Evolution

            //- Return the nucleation droplet and mass rates, and cirtical
            //- droplet diameter
            PtrList<volScalarField>& getNucFields();


        //  Access


        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
