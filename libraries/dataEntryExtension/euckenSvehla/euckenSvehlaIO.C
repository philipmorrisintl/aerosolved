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

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const euckenSvehla& euck
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const DataEntry<scalar>& >(euck)
            << token::SPACE << euck.coeffs_;
    }
    else
    {
        os  << static_cast<const DataEntry<scalar>& >(euck);
        os.write
        (
            reinterpret_cast<const char*>(&euck.coeffs_),
            sizeof(euck.coeffs_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const euckenSvehla&)"
    );

    return os;
}


void Foam::euckenSvehla::writeData(Ostream& os) const
{
    DataEntry<scalar>::writeData(os);

    os << token::SPACE << "(";

    forAll(coeffs_, i)
    {
        os << coeffs_[i] << token::SPACE;
    }

    os  << ")" << token::END_STATEMENT << nl;
}

// ************************************************************************* //
