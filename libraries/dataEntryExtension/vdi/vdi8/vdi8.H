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

\file vdi8.H
\brief VDI8 function

*/

#ifndef vdi8_H
#define vdi8_H

#include "DataEntry.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class vdi8;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const vdi8&
);

/*---------------------------------------------------------------------------*\
                           Class vdi8 Declaration
\*---------------------------------------------------------------------------*/

class vdi8
:
    public scalarDataEntry
{
    // Private data

        //- vdi8 coefficients - prefactor, factor and exponent
        List<scalar> coeffs_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const vdi8&);


public:

    //- Runtime type information
    TypeName("vdi8");


    // Constructors

        vdi8(const word& entryName, const dictionary& dict);

        //- Construct from components
        vdi8(const word& entryName, const List<scalar>&);

        //- Copy constructor
        vdi8(const vdi8& vdi);

        //- Construct and return a clone
        virtual tmp<scalarDataEntry> clone() const
        {
            return tmp<scalarDataEntry>(new vdi8(*this));
        }


    //- Destructor
    virtual ~vdi8();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return vdi8 value
            scalar value(const scalar x) const;

            //- Integrate between two (scalar) values
            scalar integrate(const scalar x1, const scalar x2) const;

            //- Return dimensioned constant value
            dimensioned<scalar> dimValue(const scalar) const;

            //- Integrate between two values and return dimensioned type
            dimensioned<scalar> dimIntegrate
            (
                const scalar x1,
                const scalar x2
            ) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const vdi8& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
