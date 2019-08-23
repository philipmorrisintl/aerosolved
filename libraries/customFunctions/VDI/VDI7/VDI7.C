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

#include "VDI7.H"
#include "mathematicalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI7<Type>::VDI7
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
            << "The VDI7 function requires 6 coefficients" << nl
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::VDI7<Type>::VDI7
(
    const VDI7<Type>& f
)
:
    Function1<Type>(f),
    c_(f.c_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI7<Type>::~VDI7()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::VDI7<Type>::value(const scalar t) const
{
    const scalar Tlim(min(t,c_[5]));
    const scalar Tr(max(1.0-Tlim/c_[5],0.0));

    return
        c_[4]
      * Foam::exp
        (
            c_[5]/Tlim
          * (
                c_[0]*Tr
              + c_[1]*Foam::pow(Tr,1.5)
              + c_[2]*Foam::pow(Tr,2.5)
              + c_[3]*Foam::pow(Tr,5.0)
            )
        );
}


template<class Type>
void Foam::Function1Types::VDI7<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << c_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
