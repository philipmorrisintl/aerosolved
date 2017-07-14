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

\file multiSpeciesCNT.H
\brief Multi-species classical nucleation theory model

A multi-species nucleation model based on the work of \cite Winkelmann:2017aa.

*/

#ifndef multiSpeciesCNT_H
#define multiSpeciesCNT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                     Class multiSpeciesCNT Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesCNT
:
    public nucleationModel
{
    // Private data

        //- Geometric standard deviation of the particle size distribution
        scalar W_;

        //- Tolerance for detecting supersaturation
        scalar TOL_;

    // Private Member Functions

        //- Disallow copy construct
        multiSpeciesCNT(const multiSpeciesCNT&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesCNT&);


protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("multiSpeciesCNT");


    // Constructors

        //- Construct from components
        multiSpeciesCNT(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~multiSpeciesCNT();


    // Member Functions

        // Evolution

            //- Return the nucleation droplet and mass rates, and cirtical
            //- droplet diameter
            PtrList<volScalarField>& getNucFields();


        //  Access

            inline const scalar& TOL() const;


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
