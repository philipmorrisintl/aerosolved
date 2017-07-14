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

/**

\file DataEntryDiffusivity.H
\brief Provides an interface to specify the diffusivity using a DataEntry
function

This method allows the binary diffusivity for a pair of selected species to be
specified using a DataEntry object. For example:

\code{.unparsed}
diffusivity
{
    Water  constant 1E-5;
}
\endcode

could be used inside the 'Air' vaporProperties part of the fluidThermoProperties
dictionary to specify a constant diffusivity of water vapor in air equal to
1E-5.

*/

#ifndef DataEntryDiffusivity_H
#define DataEntryDiffusivity_H

#include "diffusivityModel.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

// Forward declaration of classes
class DataEntryDiffusivity;

/*---------------------------------------------------------------------------*\
                    Class DataEntryDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class DataEntryDiffusivity
:
    public diffusivityModel
{
    // Private data

        //- DataEntry object
        autoPtr<DataEntry<scalar> > diffusivity_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const DataEntryDiffusivity&);


public:

    //- Runtime type information
    TypeName("DataEntryDiffusivity");


    // Constructors

        DataEntryDiffusivity
        (
            const word& entryName,
            const dictionary& dict,
            const dictionary species,
            const label a,
            const label b
        );

        //- Copy constructor
        DataEntryDiffusivity(const DataEntryDiffusivity& fsg);

        //- Construct and return a clone
        virtual tmp<diffusivityModel> clone() const
        {
            return tmp<diffusivityModel>(new DataEntryDiffusivity(*this));
        }


    //- Destructor
    virtual ~DataEntryDiffusivity();


    // Member Functions

        // Manipulation


        // Evaluation

            //- Return DataEntryDiffusivity value
            scalar value(const scalar T, const scalar p) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
