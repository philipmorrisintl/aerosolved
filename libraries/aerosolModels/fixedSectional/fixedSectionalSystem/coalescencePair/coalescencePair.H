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

#ifndef coalescencePair_H
#define coalescencePair_H

#include "secIntData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coalescencePair Declaration
\*---------------------------------------------------------------------------*/

class coalescencePair
{
protected:

    // Protected data

        //- First donor section
        label i_;

        //- Second donor section
        label j_;

        //- Mass of coalescence
        scalar s_;

        //- Sectional interpolation data
        secIntData idata_;


private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const coalescencePair&);


public:

    // Constructors

        //- Copy construct
        coalescencePair(const coalescencePair&);

        //- Construct from data
        coalescencePair
        (
            const label i,
            const label j,
            const scalar s,
            const secIntData& idata
        );


    //- Destructor
    virtual ~coalescencePair();


    // Member Functions

        // Access

            //- First donor section
            inline const label& i() const
            {
                return i_;
            }

            //- Second donor section
            inline const label& j() const
            {
                return j_;
            }

            //- Mass of coalescence
            inline const scalar& s() const
            {
                return s_;
            }

            //- Sectional interpolation data
            inline const secIntData& idata() const
            {
                return idata_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
