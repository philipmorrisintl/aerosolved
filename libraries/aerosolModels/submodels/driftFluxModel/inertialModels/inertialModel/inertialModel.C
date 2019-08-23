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

#include "inertialModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(inertialModel, 0);
defineRunTimeSelectionTable(inertialModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inertialModel::inertialModel
(
    const word& modelType,
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase(aerosol, dict, typeName, modelType),
    VMax_(dict.lookupOrDefault<scalar>("VMax", 1.0)),
    g_
    (
        IOobject
        (
            "g",
            aerosol.mesh().time().constant(),
            aerosol.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimVelocity/dimTime, vector::zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inertialModel::~inertialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField>
inertialModel::Re(const volScalarField& d, const volVectorField& V) const
{
    tmp<volScalarField> tRe
    (
        new volScalarField
        (
            IOobject
            (
                "Re",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aerosol_.mesh(),
            dimensionedScalar("Re", dimless, 0.0)
        )
    );

    const volScalarField& rhog = aerosol_.thermo().thermoCont().rho();
    const volScalarField& mug = aerosol_.thermo().thermoCont().mu();

    tRe.ref() = d*rhog*mag(V)/mug;

    return tRe;
}

tmp<volVectorField> inertialModel::limit(volVectorField& V) const
{
    tmp<volVectorField> tVs
    (
        new volVectorField
        (
            IOobject
            (
                "limit(" + V.name() + ")",
                V.mesh().time().timeName(),
                V.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            V*1.0
        )
    );

    volVectorField& Vs = tVs.ref();

    const dimensionedScalar VMin("VMin", Vs.dimensions(), SMALL);
    const dimensionedScalar VMax("VMax", Vs.dimensions(), VMax_);

    Vs = Vs/max(mag(Vs),VMin) * min(mag(Vs),VMax);

    Vs.correctBoundaryConditions();

    return tVs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
