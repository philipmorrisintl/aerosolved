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

#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo>
Foam::fluidThermo::New
(
    const fvMesh& mesh
)
{
    const word fluidThermoModel
    (
        IOdictionary
        (
            IOobject
            (
                dictName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("fluidThermoModel")
    );

    Info<< "Selecting fluid thermo model " << fluidThermoModel << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluidThermoModel);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidThermo::New"
        )   << "Unknown fluid thermo model "
            << fluidThermoModel << endl << endl
            << "Valid fluid thermo models are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<fluidThermo>(cstrIter()(mesh));
}

// ************************************************************************* //
