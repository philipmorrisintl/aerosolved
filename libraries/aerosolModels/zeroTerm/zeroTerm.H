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

#ifndef zeroTerm_H
#define zeroTerm_H

#include "aerosolModel.H"
#include "DataEntry.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{

/*---------------------------------------------------------------------------*\
                            Class zeroTerm Declaration
\*---------------------------------------------------------------------------*/

class zeroTerm
:
    public aerosolModel
{
    // Private data

        const dimensionSet dropletSizeDimension_;


    // Private Member Functions

        //- Disallow copy construct
        zeroTerm(const zeroTerm&);

        //- Disallow default bitwise assignment
        void operator=(const zeroTerm&);


protected:

    // Protected data

        const fvMesh& mesh_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("zeroTerm");


    // Constructors

        //- Construct from components
        zeroTerm(const fvMesh& mesh);


    //- Destructor
    virtual ~zeroTerm();


    // Member Functions

        // Access

            //- This is a moment method
            label modType() const {return MOMENTAEROSOLMODEL;};

            //- Return the size distribution type
            label sizeDistType() const {return NOSIZEDIST;};

            //- Droplet size dimension
            const dimensionSet& dropletSizeDimension() const
            {
                return dropletSizeDimension_;
            };

        // Functions

            void update();

            void fractionalStepInternal();
            void fractionalStepExternal();

            void checkConsistency();

            void correctSizeDistribution();

            tmp<volScalarField> dcm();

            tmp<volScalarField> dmm();

        // I-O

            //- Read info from dictionaries
            bool read();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aerosolModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zeroTermI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
