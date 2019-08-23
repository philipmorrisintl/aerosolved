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

#ifndef conData_H
#define conData_H

#include "scalarList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class conData;

// Friend operators

conData operator*(const scalar&, const conData&);
conData operator*(const conData&, const scalar&);

/*---------------------------------------------------------------------------*\
                       Class conData Declaration
\*---------------------------------------------------------------------------*/

class conData
{
protected:

    // Protected data

        //- Source coefficients
        scalarList source_;

        //- Sink coefficients
        scalarList sink_;

        //- Active
        Switch active_;


public:

    // Constructors

        //- Copy construct
        conData(const conData&);

        //- Construct from number of species
        conData(const label N);

    //- Destructor
    virtual ~conData();


    // Member Functions

        // Access

            //- Source coefficients
            inline scalarList& source()
            {
                return source_;
            }

            inline const scalarList& source() const
            {
                return source_;
            }

            //- Sink coefficients
            inline scalarList& sink()
            {
                return sink_;
            }

            inline const scalarList& sink() const
            {
                return sink_;
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


    // Member Operators

        void operator*=(const scalar&);
        void operator/=(const scalar&);
        void operator=(const conData&);


    // Friend Operators

        friend conData operator*
        (
            const scalar&,
            const conData&
        );

        friend conData operator*
        (
            const conData&,
            const scalar&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
