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

#include "exponentialLog.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialLog, 0);
    addToRunTimeSelectionTable(scalarDataEntry, exponentialLog, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponentialLog::exponentialLog(const word& entryName, const dictionary& dict)
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
            "Foam::exponentialLog::exponentialLog(const word&, const dictionary&)"
        )   << "exponentialLog coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::exponentialLog::exponentialLog
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
            "Foam::exponentialLog::exponentialLog"
            "(const word&, const List<Tuple2<scalar, scalar> >&)"
        )   << "exponentialLog coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


Foam::exponentialLog::exponentialLog(const exponentialLog& expoLog)
:
    scalarDataEntry(expoLog),
    coeffs_(expoLog.coeffs_),
    dimensions_(expoLog.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::exponentialLog::~exponentialLog()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::exponentialLog::convertTimeBase(const Time& t)
{
    scalar value = coeffs_[0];
    coeffs_[0] = t.userTimeToTime(value);
}


Foam::scalar Foam::exponentialLog::value(const scalar x) const
{
    return coeffs_[0] *
           exp
           (
               coeffs_[1] +
               coeffs_[2] * pow(x, coeffs_[3]) +
               coeffs_[4] * pow(log(x), coeffs_[5])
           );
}


Foam::scalar Foam::exponentialLog::integrate(const scalar x1, const scalar x2) const
{
    // placeholder for now, quite complex integral
    return 0.0;
}


Foam::dimensioned<Foam::scalar> Foam::exponentialLog::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::exponentialLog::dimIntegrate
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
