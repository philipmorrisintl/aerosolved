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

#include "vdi3.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vdi3, 0);
    addToRunTimeSelectionTable(scalarDataEntry, vdi3, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vdi3::vdi3(const word& entryName, const dictionary& dict)
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
            "Foam::vdi3::vdi3(const word&, const dictionary&)"
        )   << "vdi3 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::vdi3::vdi3
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
            "Foam::vdi3::vdi3"
            "(const word&, const List< <scalar> >&)"
        )   << "vdi3 coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::vdi3::vdi3(const vdi3& vdi)
:
    scalarDataEntry(vdi),
    coeffs_(vdi.coeffs_),
    dimensions_(vdi.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vdi3::~vdi3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vdi3::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}


Foam::scalar Foam::vdi3::value(const scalar T) const
{
    const scalar& A = coeffs_[0];
    const scalar& B = coeffs_[1];
    const scalar& C = coeffs_[2];
    const scalar& D = coeffs_[3];
    const scalar& E = coeffs_[4];

    return max
    (
        A + B*T + C*pow(T,2.) + D*pow(T,3.) + E*pow(T,4.),
        0.0
    );
}


Foam::scalar Foam::vdi3::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}


Foam::dimensioned<Foam::scalar> Foam::vdi3::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::vdi3::dimIntegrate
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
