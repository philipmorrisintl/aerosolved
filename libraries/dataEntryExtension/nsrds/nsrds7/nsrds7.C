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

#include "nsrds7.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nsrds7, 0);
    addToRunTimeSelectionTable(scalarDataEntry, nsrds7, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nsrds7::nsrds7(const word& entryName, const dictionary& dict)
:
    scalarDataEntry(entryName),
    coeffs_(),
    dimensions_(dimless)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);

    token firstToken(is);
    is.putBack(firstToken);
    if (firstToken == token::BEGIN_SQR)
    {
        is  >> this->dimensions_;
    }

    is  >> coeffs_;

    if (!coeffs_.size())
    {
        FatalErrorIn
        (
            "Foam::nsrds7::nsrds7(const word&, const dictionary&)"
        )   << "nsrds7 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::nsrds7::nsrds7
(
    const word& entryName,
    const List<scalar>& coeffs
)
:
    scalarDataEntry(entryName),
    coeffs_(coeffs),
    dimensions_(dimless)
{
    if (!coeffs_.size())
    {
        FatalErrorIn
        (
            "Foam::nsrds7::nsrds7"
            "(const word&, const List< <scalar> >&)"
        )   << "nsrds7 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::nsrds7::nsrds7(const nsrds7& nsrds)
:
    scalarDataEntry(nsrds),
    coeffs_(nsrds.coeffs_),
    dimensions_(nsrds.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nsrds7::~nsrds7()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nsrds7::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}


Foam::scalar Foam::nsrds7::value(const scalar x) const
{
    //a_ + b_*sqr((c_/T)/sinh(c_/T)) + d_*sqr((e_/T)/cosh(e_/T));
  return max( coeffs_[0] + coeffs_[1]*sqr((coeffs_[2]/x)/sinh(coeffs_[2]/x)) + coeffs_[3]*sqr((coeffs_[4]/x)/cosh(coeffs_[4]/x)), 0.0);
}


Foam::scalar Foam::nsrds7::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}


Foam::dimensioned<Foam::scalar> Foam::nsrds7::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::nsrds7::dimIntegrate
(
    const scalar x1,
    const scalar x2
) const
{
    return dimensioned<scalar>
    (
        "dimensionedValue",
        dimensions_,
        integrate(x1, x2)
    );
}

// ************************************************************************* //
