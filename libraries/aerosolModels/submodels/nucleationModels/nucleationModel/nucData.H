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

#ifndef nucData_H
#define nucData_H

#include "scalarList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nucData Declaration
\*---------------------------------------------------------------------------*/

class nucData
{
protected:

    // Protected data

        //- Mass fractions in the critical cluster
        scalarList z_;

        //- Size of the critical cluster
        scalar s_;

        //- Nucleation rate
        scalar J_;

        //- Active
        Switch active_;


private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const nucData&);


public:

    // Constructors

        //- Copy construct
        nucData(const nucData&);

        //- Construct from number of species
        nucData(const label N);

    //- Destructor
    virtual ~nucData();


    // Member Functions

        // Access

            //- Mass fractions in the critical cluster
            inline scalarList& z()
            {
                return z_;
            }

            inline const scalarList& z() const
            {
                return z_;
            }

            //- Size of the critical cluster
            inline scalar& s()
            {
                return s_;
            }

            inline const scalar& s() const
            {
                return s_;
            }

            //- Nucleation rate
            inline scalar& J()
            {
                return J_;
            }

            inline const scalar& J() const
            {
                return J_;
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
