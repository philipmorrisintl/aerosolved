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

\file massFlux.H
\brief functionObject to compute the mass flux at a patch of faceZone

*/

#ifndef functionObjects_massFlux_H
#define functionObjects_massFlux_H

#include "sampleFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class massFlux Declaration
\*---------------------------------------------------------------------------*/

class massFlux
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
        massFlux(const massFlux&);

        //- Disallow default bitwise assignment
        void operator=(const massFlux&);


public:

    //- Runtime type information
    TypeName("massFlux");


    // Constructors

        //- Construct from Time and dictionary
        massFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~massFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
