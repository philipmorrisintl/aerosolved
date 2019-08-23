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

#include "section.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

section::section
(
    const label i,
    const word name,
    const scalar x,
    const scalar yl,
    const scalar yu,
    const dimensionSet sizeDims,
    const fvMesh& mesh
)
:
    sectionNum_(i),
    sectionName_(name),
    x_(x),
    yl_(yl),
    yu_(yu),
    sizeDims_(sizeDims),
    M_(),
    phiInertial_
    (
        IOobject
        (
            IOobject::groupName("phiInertial", sectionName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    phiBrownian_
    (
        IOobject
        (
            IOobject::groupName("phiBrownian", sectionName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::flux
        (
            volVectorField
            (
                IOobject
                (
                    "dummy",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "dummy",
                    dimVelocity*dimDensity,
                    vector::zero
                )
            )
        )
    ),
    D_
    (
        IOobject
        (
            IOobject::groupName("D", sectionName_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("D", dimArea/dimTime, 0.0)
    ),
    validM_(false)
{
    IOobject fieldHeader
    (
        IOobject::groupName("M", sectionName_),
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    IOobject defaultFieldHeader
    (
        "M",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (fieldHeader.typeHeaderOk<volScalarField>(true))
    {
        M_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("M", sectionName_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        validM_ = true;
    }
    else if (defaultFieldHeader.typeHeaderOk<volScalarField>(true))
    {
        tmp<volScalarField> tdefaultField
        (
            new volScalarField
            (
                IOobject
                (
                    "M",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

        M_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("M", sectionName_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                tdefaultField()
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Field " << IOobject::groupName("M", sectionName_)
            << " not found"
            << " (not the field itself, nor the default field)"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

section::~section()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void section::initM(const word& patchName)
{
    M_->correctBoundaryConditions();

    label patchI = M_().mesh().boundaryMesh().findPatchID(patchName);

    if(patchI < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName << " does not exist in "
            << M_().mesh().boundaryMesh().names() << endl
            << exit(FatalError);
    }

    const fvPatchField<scalar>& patchField = M_().boundaryField()[patchI];

    scalar avgValue =
        gSum(patchField*patchField.patch().magSf())
      / gSum(patchField.patch().magSf());

    Info<< "Setting " << M_->name() << " from value on patch "
        << patchName << " to " << avgValue << endl;

    M().primitiveFieldRef() = avgValue;

    M() == dimensionedScalar("avg",M().dimensions(),avgValue);

    validM_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
