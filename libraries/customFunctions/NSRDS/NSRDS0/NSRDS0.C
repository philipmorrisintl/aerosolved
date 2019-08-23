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

#include "NSRDS0.H"
#include "mathematicalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDS0<Type>::NSRDS0
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);
    is >> c_;

    if (c_.size() != 6)
    {
        FatalErrorInFunction
            << "The NSRDS0 function requires 6 coefficients" << nl
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::NSRDS0<Type>::NSRDS0
(
    const NSRDS0<Type>& f
)
:
    Function1<Type>(f),
    c_(f.c_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDS0<Type>::~NSRDS0()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::NSRDS0<Type>::value(const scalar t) const
{
    return max
    (
        (
            (
                (
                    (
                        c_[5]
                      * t
                      + c_[4]
                    )
                  * t
                  + c_[3]
                )
              * t
              + c_[2]
            )
          * t
          + c_[1]
        )
      * t
      + c_[0],
        0.0
    );
}


template<class Type>
void Foam::Function1Types::NSRDS0<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << c_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
