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

#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diffusivityModel> Foam::diffusivityModel::New
(
    const word& entryName,
    const dictionary& dict,
    aerosolThermo& thermo,
    const label j,
    const label k
)
{
    Istream& is(dict.lookup(entryName, false));

    token firstToken(is);

    word diffusivityModelType;

    diffusivityModelType = firstToken.wordToken();

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diffusivityModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        is.putBack(firstToken);

        cstrIter = dictionaryConstructorTablePtr_->find("function");
    }

    return cstrIter()(entryName, dict, thermo, j, k);
}


// ************************************************************************* //
