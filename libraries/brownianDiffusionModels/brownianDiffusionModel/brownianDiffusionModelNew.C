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

#include "brownianDiffusionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::brownianDiffusionModel>
Foam::brownianDiffusionModel::New
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
{
    const word brownianDiffusionModelName
    (
        IOdictionary
        (
            IOobject
            (
                "brownianDiffusionProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("brownianDiffusionModel")
    );

    Info<< "Selecting Brownian diffusion model " << brownianDiffusionModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(brownianDiffusionModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "brownianDiffusionModel::New"
        )   << "Unknown diffusivity model "
            << brownianDiffusionModelName << endl << endl
            << "Valid diffusivity models are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<brownianDiffusionModel>(cstrIter()
        (mesh, aerosol));
}

// ************************************************************************* //
