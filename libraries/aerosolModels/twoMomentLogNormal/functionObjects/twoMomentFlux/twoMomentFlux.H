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

/**

\file twoMomentFlux.H
\brief functionObject to determine the number concentration flux across a patch
or faceZone

*/

#ifndef functionObjects_twoMomentFlux_H
#define functionObjects_twoMomentFlux_H

#include "sampleFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class twoMomentFlux Declaration
\*---------------------------------------------------------------------------*/

class twoMomentFlux
:
    public sampleFlux
{
protected:

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        twoMomentFlux(const twoMomentFlux&);

        //- Disallow default bitwise assignment
        void operator=(const twoMomentFlux&);


public:

    //- Runtime type information
    TypeName("twoMomentFlux");


    // Constructors

        //- Construct from Time and dictionary
        twoMomentFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~twoMomentFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
