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

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::viscosityModel>
Foam::viscosityModel::New
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
{
    const word viscosityModelName
    (
        IOdictionary
        (
            IOobject
            (
                "viscosityProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("viscosityModel")
    );

    Info<< "Selecting viscosity model " << viscosityModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(viscosityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "viscosityModel::New"
        )   << "Unknown viscosity model "
            << viscosityModelName << endl << endl
            << "Valid viscosity models are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<viscosityModel>(cstrIter()
        (mesh, thermo));
}

// ************************************************************************* //
