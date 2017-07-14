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

\file LeeChen.H
\brief Coalescence rate based on the Lee-Chen model

See \cite Lee:1984aa. This model returns the coalescence rate based on the model
given in \cite Lee:1984aa.

*/

#ifndef LeeChen_H
#define LeeChen_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                          Class LeeChen Declaration
\*---------------------------------------------------------------------------*/

class LeeChen
:
    public coalescenceModel
{
    // Private data

        const volScalarField* muEffPtr_;

    // Private Member Functions

        //- Disallow copy construct
        LeeChen(const LeeChen&);

        //- Disallow default bitwise assignment
        void operator=(const LeeChen&);


protected:

    // Protected data
        //- Geometric standard deviation of the particle size distribution
        scalar W_;

public:

    //- Runtime type information
    TypeName("LeeChen");


    // Constructors

        //- Construct from components
        LeeChen(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~LeeChen();


    // Member Functions

        // Evolution

            scalar getCoaRateCell
            (
                const scalar,
                const scalar,
                const label,
                const scalar
            );

        // Access

            inline const scalar& W() const;

        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
