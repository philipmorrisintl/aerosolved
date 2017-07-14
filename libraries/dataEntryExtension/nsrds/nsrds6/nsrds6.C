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

#include "nsrds6.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nsrds6, 0);
    addToRunTimeSelectionTable(scalarDataEntry, nsrds6, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nsrds6::nsrds6(const word& entryName, const dictionary& dict)
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
            "Foam::nsrds6::nsrds6(const word&, const dictionary&)"
        )   << "nsrds6 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::nsrds6::nsrds6
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
            "Foam::nsrds6::nsrds6"
            "(const word&, const List< <scalar> >&)"
        )   << "nsrds6 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::nsrds6::nsrds6(const nsrds6& nsrds)
:
    scalarDataEntry(nsrds),
    coeffs_(nsrds.coeffs_),
    dimensions_(nsrds.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nsrds6::~nsrds6()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nsrds6::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}


Foam::scalar Foam::nsrds6::value(const scalar x) const
{
   //Tr = T/Tc_;
   //a_*pow(1 - Tr, ((e_*Tr + d_)*Tr + c_)*Tr + b_);

    scalar Tr = min(x/coeffs_[5], 1.0);
    scalar val = 0;
    val =  max( coeffs_[0]*pow( 1 - Tr,
                (( coeffs_[4]*Tr +
                   coeffs_[3])*Tr + coeffs_[2])*Tr +
                coeffs_[1]), 0.0);
    return val;
}


Foam::scalar Foam::nsrds6::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}


Foam::dimensioned<Foam::scalar> Foam::nsrds6::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::nsrds6::dimIntegrate
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
