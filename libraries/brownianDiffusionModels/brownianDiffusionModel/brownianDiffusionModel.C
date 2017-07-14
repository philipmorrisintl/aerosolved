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

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(brownianDiffusionModel, 0);
    defineRunTimeSelectionTable(brownianDiffusionModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::brownianDiffusionModel::brownianDiffusionModel
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    IOdictionary
    (
        IOobject
        (
            "brownianDiffusionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    aerosol_(aerosol),
    thermo_(aerosol.thermo()),
    coeffs_(subDict("Coeffs")),
    DM_(0)
{
    if (aerosol_.modType() == MOMENTAEROSOLMODEL)
    {
        DM_.setSize(1);

        DM_.set
        (
            0,
            new volScalarField
            (
                IOobject
                (
                    "DM",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DM", dimLength*dimLength/dimTime, 0.0)
            )
        );
    }
    else if (aerosol_.modType() == SECTIONALAEROSOLMODEL)
    {
        DM_.setSize(aerosol_.P());

        forAll(aerosol_.x(), i)
        {
            Foam::string name("D" + aerosol_.M()[i].name());

            DM_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("DM", dimLength*dimLength/dimTime, 0.0)
                )
            );
        }
    }
    else
    {
        FatalErrorIn("Foam::brownianDiffusionModel::brownianDiffusionModel(...)")
            << "Invalid aerosol model type" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::brownianDiffusionModel::~brownianDiffusionModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
