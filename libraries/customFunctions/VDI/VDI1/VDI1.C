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

#include "VDI1.H"
#include "mathematicalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI1<Type>::VDI1
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
            << "The VDI1 function requires 6 coefficients" << nl
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::VDI1<Type>::VDI1
(
    const VDI1<Type>& f
)
:
    Function1<Type>(f),
    c_(f.c_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI1<Type>::~VDI1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::VDI1<Type>::value(const scalar t) const
{
    const scalar Tr(min(max(1.0-t/c_[5],0.0),1.0));

    return
        c_[0]*Foam::pow(Tr,0.35)
      + c_[1]*Foam::pow(Tr,2.0/3.0)
      + c_[2]*Tr
      + c_[3]*Foam::pow(Tr,4.0/3.0)
      + c_[4];
}


template<class Type>
void Foam::Function1Types::VDI1<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << c_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
