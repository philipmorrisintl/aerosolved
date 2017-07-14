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

#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diffusivityModel> Foam::diffusivityModel::New
(
    const word& entryName,
    const dictionary& dict,
    const dictionary species,
    const label a,
    const label b
)
{
    Istream& is(dict.lookup(entryName, false));

    token firstToken(is);

    word diffusivityModelType = firstToken.wordToken();

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diffusivityModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        // If this diffusivity model doesn't exist, treat as DataEntry

        cstrIter = dictionaryConstructorTablePtr_->find("DataEntryDiffusivity");
    }

    return autoPtr<diffusivityModel>(cstrIter()(entryName, dict, species, a, b));
}


// ************************************************************************* //
