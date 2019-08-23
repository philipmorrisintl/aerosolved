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

\file sectionalFlux.H
\brief functionObject to determine the sectional flux at a patch or faceZone

*/

#ifndef functionObjects_sectionalFlux_H
#define functionObjects_sectionalFlux_H

#include "sampleFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class sectionalFlux Declaration
\*---------------------------------------------------------------------------*/

class sectionalFlux
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
        sectionalFlux(const sectionalFlux&);

        //- Disallow default bitwise assignment
        void operator=(const sectionalFlux&);


public:

    //- Runtime type information
    TypeName("sectionalFlux");


    // Constructors

        //- Construct from Time and dictionary
        sectionalFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~sectionalFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
