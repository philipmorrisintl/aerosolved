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

\file sutherland.H
\brief Sutherland function

*/

#ifndef sutherland_H
#define sutherland_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sutherland;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const sutherland&
);

/*---------------------------------------------------------------------------*\
                        Class sutherland Declaration
\*---------------------------------------------------------------------------*/

class sutherland
:
    public scalarDataEntry
{
    // Private data

        //- sutherland coefficients - prefactor, factor and suthernent
        List<scalar> coeffs_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const sutherland&);


public:

    //- Runtime type information
    TypeName("sutherland");


    // Constructors

        sutherland(const word& entryName, const dictionary& dict);

        //- Construct from components
        sutherland(const word& entryName, const List<scalar>&);

        //- Copy constructor
        sutherland(const sutherland& suther);

        //- Construct and return a clone
        virtual tmp<scalarDataEntry> clone() const
        {
            return tmp<scalarDataEntry>(new sutherland(*this));
        }


    //- Destructor
    virtual ~sutherland();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return sutherland value
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
            const sutherland& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
