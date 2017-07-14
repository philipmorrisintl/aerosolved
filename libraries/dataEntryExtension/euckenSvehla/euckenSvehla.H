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

\file euckenSvehla.H
\brief Eucken-Svehla function

*/

#ifndef euckenSvehla_H
#define euckenSvehla_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class euckenSvehla;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const euckenSvehla&
);

/*---------------------------------------------------------------------------*\
                        Class euckenSvehla Declaration
\*---------------------------------------------------------------------------*/

class euckenSvehla
:
    public scalarDataEntry
{
    // Private data

        //- euckenSvehla coefficients - prefactor, factor and eucknent
        List<scalar> coeffs_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const euckenSvehla&);

        //- Compute janaf step value
        scalar euckenSvehlaValue(const scalar x, label i) const;



public:

    //- Runtime type information
    TypeName("euckenSvehla");


    // Constructors

        euckenSvehla(const word& entryName, const dictionary& dict);

        //- Construct from components
        euckenSvehla(const word& entryName, const List<scalar>&);

        //- Copy constructor
        euckenSvehla(const euckenSvehla& euck);

        //- Construct and return a clone
        virtual tmp<scalarDataEntry> clone() const
        {
            return tmp<scalarDataEntry>(new euckenSvehla(*this));
        }


    //- Destructor
    virtual ~euckenSvehla();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return euckenSvehla value
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
            const euckenSvehla& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
