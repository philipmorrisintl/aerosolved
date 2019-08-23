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

#include "VDI10.H"
#include "mathematicalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI10<Type>::VDI10
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

    if (c_.size() != 8)
    {
        FatalErrorInFunction
            << "The VDI10 function requires 8 coefficients" << nl
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::VDI10<Type>::VDI10
(
    const VDI10<Type>& f
)
:
    Function1<Type>(f),
    c_(f.c_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::VDI10<Type>::~VDI10()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::VDI10<Type>::value(const scalar t) const
{
    const scalar R(constant::physicoChemical::R.value());
    const scalar Tr(t/(c_[0]+t));

    return max
    (
        R/(c_[7]*1e-3)
      * (
            c_[1]
          + (
                (c_[2]-c_[1])
              * Foam::pow(Tr,2.0)
            )
          * (
                1.0
              - c_[0]/(c_[0]+t)
              * (
                    c_[3]
                  + c_[4]*Tr
                  + c_[5]*pow(Tr,2.0)
                  + c_[6]*pow(Tr,3.0)
                )
            )
        ),
        0.0
    );
}


template<class Type>
void Foam::Function1Types::VDI10<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << c_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
