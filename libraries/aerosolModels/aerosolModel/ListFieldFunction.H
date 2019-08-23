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

#ifndef ListFieldFunction_H
#define ListFieldFunction_H

#include "PtrList.H"
#include "List.H"
#include "Field.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
List<Type> entryList
(
    const PtrList<Field<Type>>& f,
    const label& i
)
{
    List<Type> s(f.size());

    forAll(f, j)
    {
        s[j] = f[j][i];
    }

    return s;
}

template<class Type>
List<Type> entryList
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh>>& f,
    const label& i
)
{
    List<Type> s(f.size());

    forAll(f, j)
    {
        s[j] = f[j].field()[i];
    }

    return s;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
