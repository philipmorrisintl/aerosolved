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

#include "condensationModel.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(condensationModel, 0);
defineRunTimeSelectionTable(condensationModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

condensationModel::condensationModel
(
    const word& modelType,
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    aerosolSubModelBase(aerosol, dict, typeName, modelType),
    activity_()
{
    if (modelType != "none")
    {
        activity_ =
            activityCoeffModel::New(aerosol, dict.subDict("activityCoeff"));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

condensationModel::~condensationModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> condensationModel::Qdot
(
    const PtrList<volScalarField>& I
) const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                aerosol_.mesh().time().timeName(),
                aerosol_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aerosol_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if
    (
        modelType() != "none"
     && dict().subDict("heatOfVaporization").lookup("active")
    )
    {
        scalarField& Qdot = tQdot.ref();

        const rhoAerosolPhaseThermo& thermoDisp = aerosol_.thermo().thermoDisp();

        const speciesTable& activeSpecies = aerosol_.thermo().activeSpecies();

        forAll(activeSpecies, j)
        {
            const scalarField& Ij = I[j].field();

            const scalar hj = thermoDisp.composition().Hc(j);

            forAll(Qdot, celli)
            {
                Qdot[celli] += hj*Ij[celli];
            }
        }

        const dictionary& hVapDict = dict().subDict("heatOfVaporization");

        if (hVapDict.lookupOrDefault("limit", false))
        {
            const scalarField& rDeltaT = aerosol_.getRDeltaT();
            const scalarField Cv(aerosol_.thermo().Cv());
            const scalarField& rho = aerosol_.rho().field();

            const scalar maxDeltaTemp
            (
                readScalar
                (
                    dict().subDict("heatOfVaporization").lookup("maxDeltaTemp")
                )
            );

            forAll(Qdot, celli)
            {
                Qdot[celli] =
                    sign(Qdot[celli])
                  * min
                    (
                        mag(Qdot[celli]),
                        maxDeltaTemp*rDeltaT[celli]*Cv[celli]*rho[celli]
                    );
            }
        }
    }

    return tQdot;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
