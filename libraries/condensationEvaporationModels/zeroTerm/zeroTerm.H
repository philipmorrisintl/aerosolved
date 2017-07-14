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

#ifndef zeroTerm_H
#define zeroTerm_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{

/*---------------------------------------------------------------------------*\
                            Class zeroTerm Declaration
\*---------------------------------------------------------------------------*/

class zeroTerm
:
    public condensationEvaporationModel
{
    // Private data

    // Private Member Functions

        //- Disallow copy construct
        zeroTerm(const zeroTerm&);

        //- Disallow default bitwise assignment
        void operator=(const zeroTerm&);


protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("zeroTerm");


    // Constructors

        //- Construct from components
        zeroTerm(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~zeroTerm();


    // Member Functions

        // Evolution

            List< List<scalar> > getCondRateListCell
            (
                const List<scalar>& z,
                const label jCell
            );

            List< List<scalar> > getEtaGammaListCell
            (
                const List<scalar>& z,
                const label jCell
            );

            //- Return droplet diameter for a given list of zetas and cell index
            List<scalar> psiInv
            (
                const List<scalar>& zeta,
                const label jCell
            );

        //  Access

        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace condensationEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
