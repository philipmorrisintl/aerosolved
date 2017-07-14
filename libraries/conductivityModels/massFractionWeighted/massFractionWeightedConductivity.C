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

#include "massFractionWeightedConductivity.H"
#include "makeConductivityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace conductivityModels
{
    makeconductivityModelTypes(massFractionWeightedConductivity, conductivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityModels::massFractionWeightedConductivity::massFractionWeightedConductivity
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    conductivityModel(mesh, thermo)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityModels::massFractionWeightedConductivity::~massFractionWeightedConductivity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conductivityModels::massFractionWeightedConductivity::update()
{
    // Determine kEff by using temperature and mass fractions

    scalarField& kEffCells = kEff_.internalField();
    const scalarField& TCells = thermo_.T().internalField();

    PtrList<DataEntry<scalar> > dataEntriesListK_v =
        thermo_.getProperty("k", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> > dataEntriesListK_l =
        thermo_.getProperty("k", fluidThermo::LIQUID);

    kEffCells = 0.0;

    forAll(thermo_.species(), i)
    {
        const volScalarField& Yi = thermo_.Y()[i];
        const volScalarField& Zi = thermo_.Z()[i];

        forAll(kEffCells, cellj)
        {
            kEffCells[cellj] +=
                max(Yi.internalField()[cellj],0.0) * dataEntriesListK_v[i].value(TCells[cellj])
              + max(Zi.internalField()[cellj],0.0) * dataEntriesListK_l[i].value(TCells[cellj]);
        }
    }

    // Boundary fields

    forAll(kEff_.boundaryField(), patchi)
    {
        const scalarField& T = thermo_.T().boundaryField()[patchi];
        scalarField& kEff = kEff_.boundaryField()[patchi];

        kEff = 0.0;

        forAll(thermo_.species(), i)
        {
            const scalarField& Yi = thermo_.Y()[i].boundaryField()[patchi];
            const scalarField& Zi = thermo_.Z()[i].boundaryField()[patchi];

            forAll(kEff, facej)
            {
                kEff[facej] +=
                    max(Yi[facej],0.0) * dataEntriesListK_v[i].value(T[facej])
                  + max(Zi[facej],0.0) * dataEntriesListK_l[i].value(T[facej]);
            }
        }
    }
}

bool Foam::conductivityModels::massFractionWeightedConductivity::read()
{
    thermo_.readPropertyBoth("k", thermo_.species());

    return true;
}

// ************************************************************************* //
