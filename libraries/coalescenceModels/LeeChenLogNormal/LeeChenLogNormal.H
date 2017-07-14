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

\file LeeChenLogNormal.H
\brief Coalescence rate based on the Lee-Chen model and assuming a log-normal
distribution.

See \cite Lee:1984aa. This model returns an effective coalescence rate assuming
log-normality of the size distribution. This model checks if the aerosolModel
is of type 'MOMENTAEROSOLMODEL'.

*/

#ifndef LeeChenLogNormal_H
#define LeeChenLogNormal_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                            Class LeeChenLogNormal Declaration
\*---------------------------------------------------------------------------*/

class LeeChenLogNormal
:
    public coalescenceModel
{
    // Private data

    // Private Member Functions

        //- Disallow copy construct
        LeeChenLogNormal(const LeeChenLogNormal&);

        //- Disallow default bitwise assignment
        void operator=(const LeeChenLogNormal&);


protected:

    // Protected data
        //- Geometric standard deviation of the particle size distribution
        scalar W_;

public:

    //- Runtime type information
    TypeName("LeeChenLogNormal");


    // Constructors

        //- Construct from components
        LeeChenLogNormal(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~LeeChenLogNormal();


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
