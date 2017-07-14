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

#include "vdi9.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vdi9, 0);
    addToRunTimeSelectionTable(scalarDataEntry, vdi9, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vdi9::vdi9(const word& entryName, const dictionary& dict)
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
            "Foam::vdi9::vdi9(const word&, const dictionary&)"
        )   << "vdi9 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::vdi9::vdi9
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
            "Foam::vdi9::vdi9"
            "(const word&, const List< <scalar> >&)"
        )   << "vdi9 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::vdi9::vdi9(const vdi9& vdi)
:
    scalarDataEntry(vdi),
    coeffs_(vdi.coeffs_),
    dimensions_(vdi.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vdi9::~vdi9()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vdi9::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}

Foam::scalar Foam::vdi9::value(const scalar T) const
{
    const scalar& A = coeffs_[0];
    const scalar& B = coeffs_[1];
    const scalar& C = coeffs_[2];
    const scalar& D = coeffs_[3];
    const scalar& E = coeffs_[4];
    const scalar& M = coeffs_[5];
    const scalar& Tc = coeffs_[6];
    const scalar R = 8.314462;

    scalar Tr = max(1.0 - T/Tc, 0.0);

    return max
    (
        R/(M*1e-3)*Tc * ( A*pow(Tr,1./3.) + B*pow(Tr,2./3.) + C*Tr
      + D*pow(Tr,2.0) + E*pow(Tr,6.0) ),
        0.0
    );
}

Foam::scalar Foam::vdi9::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}

Foam::dimensioned<Foam::scalar> Foam::vdi9::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}

Foam::dimensioned<Foam::scalar> Foam::vdi9::dimIntegrate
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
