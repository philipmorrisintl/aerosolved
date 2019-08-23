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

#ifndef secIntData_H
#define secIntData_H

#include "labelList.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class secIntData Declaration
\*---------------------------------------------------------------------------*/

class secIntData
{
protected:

    // Protected data

        //- Section indices
        labelList i_;

        //- Section weights
        scalarList w_;

        //- Interpolation point
        scalar xi_;


public:

    // Constructors

        //- Construct from data
        secIntData(const labelList&, const scalarList&, const scalar&);

        //- Construct empty from number of weights
        secIntData(const label&);

    //- Destructor
    virtual ~secIntData();


    // Member Functions

        // Access

            //- Indices

            inline labelList& i()
            {
                return i_;
            }

            inline const labelList& i() const
            {
                return i_;
            }

            //- Weights

            inline scalarList& w()
            {
                return w_;
            }

            inline const scalarList& w() const
            {
                return w_;
            }

            //- Interpolation point

            inline scalar& xi()
            {
                return xi_;
            }

            inline const scalar& xi() const
            {
                return xi_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
