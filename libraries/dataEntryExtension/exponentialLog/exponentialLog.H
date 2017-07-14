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

\file exponentialLog.H
\brief Exponential-log function

This DataEntry function can be used to describe the function:

\f[
    a \exp\left(b+c x^d + e \log^f(x)\right)
\f]

and can be specified in a dictionary for a variable "var" as:

\code{.unparsed}
var         exponentialLog (a b c d e f);
\endcode

*/

#ifndef exponentialLog_H
#define exponentialLog_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class exponentialLog;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const exponentialLog&
);

/*---------------------------------------------------------------------------*\
                        Class exponentialLog Declaration
\*---------------------------------------------------------------------------*/

class exponentialLog
:
    public scalarDataEntry
{
    // Private data

        //- exponentialLog coefficients - prefactor, factor and exponent
        List<scalar> coeffs_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const exponentialLog&);


public:

    //- Runtime type information
    TypeName("exponentialLog");


    // Constructors

        exponentialLog(const word& entryName, const dictionary& dict);

        //- Construct from components
        exponentialLog(const word& entryName, const List<scalar>&);

        //- Copy constructor
        exponentialLog(const exponentialLog& expoLog);

        //- Construct and return a clone
        virtual tmp<scalarDataEntry> clone() const
        {
            return tmp<scalarDataEntry>(new exponentialLog(*this));
        }


    //- Destructor
    virtual ~exponentialLog();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return exponentialLog value
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
            const exponentialLog& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
