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

#include "massFractionWeightedViscosity.H"
#include "makeViscosityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    makeViscosityModelTypes(massFractionWeightedViscosity, viscosityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::massFractionWeightedViscosity::massFractionWeightedViscosity
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    viscosityModel(mesh, thermo)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscosityModels::massFractionWeightedViscosity::~massFractionWeightedViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::viscosityModels::massFractionWeightedViscosity::update()
{
    // Determine muEff by using temperature and mass fractions

    scalarField& muEffCells = muEff_.internalField();
    const scalarField& TCells = thermo_.T().internalField();

    PtrList<DataEntry<scalar> > dataEntriesListMu_v =
        thermo_.getProperty("mu", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> > dataEntriesListMu_l =
        thermo_.getProperty("mu", fluidThermo::LIQUID);

    muEffCells = 0.0;

    forAll(thermo_.species(), i)
    {
        const volScalarField& Yi = thermo_.Y()[i];
        const volScalarField& Zi = thermo_.Z()[i];

        forAll(muEffCells, cellj)
        {
            muEffCells[cellj] +=
                max(Yi.internalField()[cellj], 0.0) * dataEntriesListMu_v[i].value(TCells[cellj])
              + max(Zi.internalField()[cellj], 0.0) * dataEntriesListMu_l[i].value(min(TCells[cellj], thermo_.Tc()[i]));
        }
    }

    // Boundary fields

    forAll(muEff_.boundaryField(), patchi)
    {
        const scalarField& T = thermo_.T().boundaryField()[patchi];
        scalarField& muEff = muEff_.boundaryField()[patchi];

        muEff = 0.0;

        forAll(thermo_.species(), i)
        {
            const scalarField& Yi = thermo_.Y()[i].boundaryField()[patchi];
            const scalarField& Zi = thermo_.Z()[i].boundaryField()[patchi];

            forAll(muEff, facej)
            {
                muEff[facej] +=
                    max(Yi[facej], 0.0) * dataEntriesListMu_v[i].value(T[facej])
                  + max(Zi[facej], 0.0) * dataEntriesListMu_l[i].value(T[facej]);
            }
        }
    }
}

bool Foam::viscosityModels::massFractionWeightedViscosity::read()
{
    thermo_.readPropertyBoth("mu");

    return true;
}


// ************************************************************************* //
