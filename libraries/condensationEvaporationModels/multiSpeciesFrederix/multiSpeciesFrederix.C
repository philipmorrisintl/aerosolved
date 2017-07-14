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

#include "condensationEvaporationModel.H"
#include "multiSpeciesFrederix.H"

#include "makeCondensationEvaporationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{
    makeCondensationEvaporationTypes(multiSpeciesFrederix, condensationEvaporationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::multiSpeciesFrederix::multiSpeciesFrederix
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    condensationEvaporationModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::multiSpeciesFrederix::~multiSpeciesFrederix()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesFrederix::getCondRateListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{

    const scalar pi = constant::mathematical::pi;

    const dictionary& speciesPhaseChange = thermo_.speciesPhaseChange();
    const label Nps = speciesPhaseChange.size();

    const dictionary& species = thermo_.species();
    const label N = species.size();

    const scalar& T = thermo_.T()[jCell];
    const scalar& rho = thermo_.rho()[jCell];
    const scalar p = thermo_.p1()[jCell] + thermo_.p0().value();

    const List<scalar>& M = thermo_.M();

    const PtrList<volScalarField>& Y = thermo_.Y();
    const PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<DataEntry<scalar> >& dataEntriesListP_s =
        thermo_.getProperty("P_s", fluidThermo::VAPOR);

    PtrList<DataEntry<scalar> >& dataEntriesListRho_l =
        thermo_.getProperty("rho", fluidThermo::LIQUID);

    // Create output list

    List<scalarList> I(z.size());

    forAll(z, i)
    {
        I[i] = scalarList(Nps, 0.0);
    }

    // Check if we have droplets

    scalar sumZ(0.0);

    forAll(speciesPhaseChange, k)
    {
        sumZ += max(Z[k][jCell],0.0);
    }

    if (sumZ > 0.0)
    {
        // Prepare some stuff

        scalarList pEq(Nps);

        scalar sumZoverRho(0.0);
        scalar sumZoverM(0.0);

        forAll(speciesPhaseChange, k)
        {
            pEq[k] =  dataEntriesListP_s[k].value(T);

            sumZoverRho += max(Z[k][jCell],0.0)
                     / stabilise(dataEntriesListRho_l[k].value(T), SMALL);
            sumZoverM += max(Z[k][jCell],0.0)/M[k];
        }

        scalar rhod = sumZ/stabilise(sumZoverRho, SMALL);

        scalar sumY(0.0);
        scalar sumYoverM(0.0);

        forAll(species, k)
        {
            sumY += max(Y[k][jCell],0.0);
            sumYoverM += max(Y[k][jCell],0.0)/M[k];
        }

        scalar Mv = sumY/sumYoverM;

        // Compute liquid mole fraction w.r.t. liquid phase

        scalarList omega(Nps);

        forAll(speciesPhaseChange, k)
        {
            omega[k] = max(Z[k][jCell],0.0)/M[k]/sumZoverM;
        }

        // Compute droplet surface equilibrium mole fraction w.r.t. gas phase
        // and the molar mean molecular weight of the equilibirum mixture

        scalarList chiEq(Nps);

        forAll(speciesPhaseChange, k)
        {
            chiEq[k] = omega[k]*pEq[k]/p;
        }

        // Compute rates

        forAll(speciesPhaseChange, j)
        {
            scalar D12 = thermo_.getDiffusivity(j, N-1, T, p);

            // Droplet surface equilibrium mass fraction w.r.t. the mixture

            scalar YEq = chiEq[j] * M[j]
                       / (chiEq[j] * M[j] + (1.0 - chiEq[j]) * Mv);

            forAll(z, i)
            {
                scalar d = pow(z[i]/rhod * 6.0/pi, 1.0/3.0);

                I[i][j] = 2.0*pi * D12 * d * rho * (max(Y[j][jCell],0.0) - YEq);
            }
        }
    }

    return I;
}

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesFrederix::getEtaGammaListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{
    const dictionary& speciesPhaseChange = thermo_.speciesPhaseChange();
    const label Nps = speciesPhaseChange.size();

    // Create output list

    List<scalarList> I(z.size());

    forAll(z, i)
    {
        I[i] = scalarList(Nps+1, 0.0);
    }

    FatalErrorIn("Foam::condensationEvaporationModels::multiSpeciesFrederix::getEtaGammaListCell(...)")
        << "Not implemented yet." << exit(FatalError);

    return I;
}

Foam::List<Foam::scalar>
Foam::condensationEvaporationModels::multiSpeciesFrederix::psiInv
(
    const Foam::List<Foam::scalar>& zeta,
    const Foam::label jCell
)
{
    List<scalar> z(zeta.size(), -1.0);

    FatalErrorIn("Foam::condensationEvaporationModels::multiSpeciesFrederix::getEtaGammaListCell(...)")
        << "Not implemented yet." << exit(FatalError);

    return z;
}

bool Foam::condensationEvaporationModels::multiSpeciesFrederix::read()
{
    if (condensationEvaporationModel::read())
    {
        thermo_.readProperty("P_s", fluidThermo::VAPOR, thermo_.speciesPhaseChange());
        thermo_.readProperty("sigma", fluidThermo::LIQUID, thermo_.speciesPhaseChange());
        thermo_.readProperty("rho", fluidThermo::LIQUID, thermo_.speciesPhaseChange());

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
