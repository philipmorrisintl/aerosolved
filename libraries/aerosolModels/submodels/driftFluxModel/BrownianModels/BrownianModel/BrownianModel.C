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

#include "BrownianModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BrownianModel, 0);
defineRunTimeSelectionTable(BrownianModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BrownianModel::BrownianModel
(
    const word& modelType,
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase(aerosol, dict, typeName, modelType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BrownianModel::~BrownianModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> BrownianModel::D(const volScalarField& d) const
{
    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedScalar("D", sqr(dimLength)/dimTime, 0.0)
        )
    );

    volScalarField& D = tD.ref();

    forAll(D.field(), celli)
    {
        D.field()[celli] = this->D(d.field()[celli], celli);
    }

    forAll(D.boundaryField(), patchi)
    {
        D.boundaryFieldRef()[patchi] =
            this->D(d.boundaryField()[patchi], patchi);
    }

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
