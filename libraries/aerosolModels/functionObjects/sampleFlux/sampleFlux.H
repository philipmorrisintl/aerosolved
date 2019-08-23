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

#ifndef functionObjects_sampleFlux_H
#define functionObjects_sampleFlux_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class sampleFlux Declaration
\*---------------------------------------------------------------------------*/

class sampleFlux
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

        //- List of patches
        wordList patchNames_;

        //- List of faceZoneSets
        wordList faceZoneSetNames_;

        //- List of patchIDs
        labelList patchIDs_;

        //- List of faceZoneSetIDs
        labelList faceZoneSetIDs_;

        //- List of flux field names
        wordList fluxFieldNames_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const = 0;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        sampleFlux(const sampleFlux&);

        //- Disallow default bitwise assignment
        void operator=(const sampleFlux&);


public:

    //- Runtime type information
    TypeName("sampleFlux");


    // Constructors

        //- Construct from Time and dictionary
        sampleFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~sampleFlux();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
