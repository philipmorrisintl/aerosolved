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

#ifndef coaData_H
#define coaData_H

#include "scalarList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coaData Declaration
\*---------------------------------------------------------------------------*/

class coaData
{
protected:

    // Protected data

        //- Polynomial weights
        scalarList w_;

        //- Polynomial p powers
        scalarList p_;

        //- Polynomial q powers
        scalarList q_;

        //- Active
        Switch active_;


private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const coaData&);


public:

    // Constructors

        //- Copy construct
        coaData(const coaData&);

        //- Construct from the number of polynomial entries
        coaData(const label N);

        //- Construct from w, p and q lists
        coaData
        (
            const scalarList& w,
            const scalarList& p,
            const scalarList& q,
            const Switch& active
        );


    //- Destructor
    virtual ~coaData();


    // Member Functions

        // Access

            //- Polynomial weights
            inline scalarList& w()
            {
                return w_;
            }

            inline const scalarList& w() const
            {
                return w_;
            }

            //- Polynomial p powers
            inline scalarList& p()
            {
                return p_;
            }

            inline const scalarList& p() const
            {
                return p_;
            }

            //- Polynomial q powers
            inline scalarList& q()
            {
                return q_;
            }

            inline const scalarList& q() const
            {
                return q_;
            }

            //- Active
            inline Switch& active()
            {
                return active_;
            }

            inline const Switch& active() const
            {
                return active_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
