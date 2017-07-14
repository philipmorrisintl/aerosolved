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

#include "euckenSvehla.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(euckenSvehla, 0);
    addToRunTimeSelectionTable(scalarDataEntry, euckenSvehla, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::euckenSvehla::euckenSvehla(const word& entryName, const dictionary& dict)
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
            "Foam::euckenSvehla::euckenSvehla(const word&, const dictionary&)"
        )   << "euckenSvehla coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::euckenSvehla::euckenSvehla
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
            "Foam::euckenSvehla::euckenSvehla"
            "(const word&, const List<Tuple2<scalar, scalar> >&)"
        )   << "euckenSvehla coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::euckenSvehla::euckenSvehla(const euckenSvehla& euck)
:
    scalarDataEntry(euck),
    coeffs_(euck.coeffs_),
    dimensions_(euck.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::euckenSvehla::~euckenSvehla()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::euckenSvehla::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}


Foam::scalar Foam::euckenSvehla::value(const scalar x) const
{
    //scalar mu = coeffs_[3] * pow(x , 3.0/2.0)/( x + coeffs_[4]); // Sutherland
    //scalar cp = euckenSvehlaValue(x,5); // janaf

    scalar mu =  coeffs_[4]*pow(x, coeffs_[5])/(1.0 + coeffs_[6]/x + coeffs_[7]/sqr(x)) ;  // nsrds2
    scalar cp =  coeffs_[8] + coeffs_[9]*sqr((coeffs_[10]/x)/sinh(coeffs_[10]/x)) + coeffs_[11]*sqr((coeffs_[12]/x)/cosh(coeffs_[12]/x)); // nsrds7

    return (1.32+1.77*coeffs_[0]*coeffs_[2]/coeffs_[1]/1.0e-3/cp) * cp * mu / coeffs_[2];
}

Foam::scalar Foam::euckenSvehla::euckenSvehlaValue(const scalar x, label i) const
{
    return (i < coeffs_.size() ? coeffs_[i] + x * euckenSvehlaValue(x, i+1) : 0.0);
}

Foam::scalar Foam::euckenSvehla::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}


Foam::dimensioned<Foam::scalar> Foam::euckenSvehla::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::euckenSvehla::dimIntegrate
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
