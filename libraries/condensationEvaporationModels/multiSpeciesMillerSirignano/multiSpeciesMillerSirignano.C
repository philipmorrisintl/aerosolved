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
#include "multiSpeciesMillerSirignano.H"

#include "makeCondensationEvaporationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{
    #ifndef DOXYGENIGNORE

    template<>
    const char* NamedEnum
    <
        SherwoodType,2
    >::names[] =
    {
        "natural",
        "forced"
    };

    template<>
    const char* NamedEnum
    <
        potentialType,3
    >::names[] =
    {
        "classical",
        "massAnalogyI",
        "massAnalogyII",
    };

    template<>
    const char* NamedEnum
    <
        diffusionType,2
    >::names[] =
    {
        "speciesSpecific",
        "mixture"
    };

    #endif

    makeCondensationEvaporationTypes(multiSpeciesMillerSirignano, condensationEvaporationModel);
}
}

const Foam::NamedEnum<Foam::condensationEvaporationModels::SherwoodType, 2>
    Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::SherwoodTypeNames;

const Foam::NamedEnum<Foam::condensationEvaporationModels::potentialType, 3>
    Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::potentialTypeNames;

const Foam::NamedEnum<Foam::condensationEvaporationModels::diffusionType, 2>
    Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::diffusionTypeNames;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::multiSpeciesMillerSirignano
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

Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::~multiSpeciesMillerSirignano()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::getCondRateListCell
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
    const scalar& mu = mesh_.lookupObject<volScalarField>("muEff")[jCell];
    const vector& U = mesh_.lookupObject<volVectorField>("U")[jCell];

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

        scalarList pSat(Nps);

        scalar sumZoverRho(0.0);
        scalar sumZoverM(0.0);
        scalar sumYps(0.0);

        forAll(speciesPhaseChange, k)
        {
            pSat[k] =  dataEntriesListP_s[k].value(T);

            sumZoverRho += max(Z[k][jCell],0.0)
                     / stabilise(dataEntriesListRho_l[k].value(T), SMALL);
            sumZoverM += max(Z[k][jCell],0.0)/M[k];

            sumYps += Y[k][jCell];
        }

        scalar rhod = sumZ/stabilise(sumZoverRho, SMALL);

        scalar sumY(0.0);
        scalar sumYoverM(0.0);

        forAll(species, k)
        {
            sumY += max(Y[k][jCell],0.0);
            sumYoverM += max(Y[k][jCell],0.0)/M[k];
        }

        scalar Mv = sumY/stabilise(sumYoverM, SMALL);

        scalarList X(N);
        scalarList W(N);

        forAll(species, k)
        {
            X[k] = max(Y[k][jCell],0.0)/M[k] / stabilise(sumYoverM, SMALL);
            W[k] = max(Z[k][jCell],0.0)/M[k] / stabilise(sumZoverM, SMALL);
        }

        // Compute liquid mole fraction w.r.t. liquid phase

        scalarList omega(Nps);

        forAll(speciesPhaseChange, k)
        {
            omega[k] = max(Z[k][jCell],0.0)/M[k]/stabilise(sumZoverM, SMALL);
        }

        // Compute droplet surface equilibrium mole fraction w.r.t. gas phase
        // and the molar mean molecular weight of the equilibirum mixture

        scalarList chiSat(Nps);

        forAll(speciesPhaseChange, k)
        {
            chiSat[k] = omega[k]*pSat[k]/p;
        }

        // Compute saturation mass fractions

        scalarList Ysat(Nps, 0.0);
        scalar sumYsat(0.0);

        forAll(speciesPhaseChange, j)
        {
            Ysat[j] = chiSat[j] * M[j]
                    / (chiSat[j] * M[j] + (1.0 - chiSat[j]) * Mv);

            sumYsat += Ysat[j];
        }

        // Droplet diameters

        scalarList d(z.size(), 0.0);

        forAll(z, i)
        {
            d[i] = pow(z[i]/rhod * 6.0/pi, 1.0/3.0);
        }

        // Spalding number

        const scalar Bj = (sumYsat - sumYps) / (1-sumYsat);

        // Species mixture diffusivities

        scalarList D(Nps, 0.0);

        forAll(speciesPhaseChange, j)
        {
            forAll(species, k)
            {
                if (k != j)
                {
                    D[j] += X[k] / thermo_.getDiffusivity(j, k, T, p);
                }
            }

            D[j] = (1.0-Y[j][jCell]) / D[j];
        }

        // Mixture mean diffusivity

        scalar DD(0.0);
        scalar sumXWps(0.0);

        forAll(speciesPhaseChange, j)
        {
            if (Bj < 0.0)
            {
                DD += (X[j]+SMALL) * D[j];
                sumXWps += (X[j]+SMALL);
            }
            else
            {
                DD += (W[j]+SMALL) * D[j];
                sumXWps += (W[j]+SMALL);
            }
        }

        DD /= sumXWps;

        // Driving potential

        scalar H(0.0);

        switch (potential_)
        {
            case CLASSICAL:
            {
                H = log(1.0+Bj);
            }
            break;

            case MASSANALOGYI:
            {
                H = Bj;
            }
            break;

            case MASSANALOGYII:
            {
                H = Bj*(1.0-sumYsat);
            }
        }

        // Sherwood numbers (size dependent)

        scalarList Sh(z.size(), 2.0);

        switch(Sherwood_)
        {
            case NATURALCONVECTION:
            {
                // Sherwood is 2 (already initialized as such)
            }
            break;

            case FORCEDCONVECTION:
            {
                const scalar Sc(mu / (rho*DD));

                forAll(z, i)
                {
                    vector V(U);

                    if (aerosol_.doDrift())
                    {
                        V = aerosol_.V()[i][jCell];
                    }

                    const scalar Re(d[i] * rho * mag(U-V) / mu);

                    const scalar F(pow(1.0+Bj, 0.7) * log(1.0+Bj) / Bj);

                    Sh[i] += 0.6*pow(Sc, 1.0/3.0) * sqrt(Re) / F;
                }
            }
            break;
        }

        // Condensation rates

        forAll(speciesPhaseChange, j)
        {
            scalar epsj(((1.0+Bj)*Ysat[j] - max(Y[j][jCell],0.0)) / Bj);
            scalar Dj(0.0);

            switch (diffusion_)
            {
                case SPECIESSPECIFIC:
                {
                    Dj = D[j];
                }
                break;

                case MIXTURE:
                {
                    Dj = DD;
                }
                break;
            }

            forAll(z, i)
            {
                I[i][j] = - pi * epsj * Sh[i] * rho * Dj * d[i] * H;
            }
        }
    }

    return I;
}

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::getEtaGammaListCell
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

    FatalErrorIn("Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::getEtaGammaListCell(...)")
        << "Not implemented yet." << exit(FatalError);

    return I;
}

Foam::List<Foam::scalar>
Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::psiInv
(
    const Foam::List<Foam::scalar>& zeta,
    const Foam::label jCell
)
{
    List<scalar> z(zeta.size(), -1.0);

    FatalErrorIn("Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::getEtaGammaListCell(...)")
        << "Not implemented yet." << exit(FatalError);

    return z;
}

bool Foam::condensationEvaporationModels::multiSpeciesMillerSirignano::read()
{
    if (condensationEvaporationModel::read())
    {
        thermo_.readProperty("P_s", fluidThermo::VAPOR, thermo_.speciesPhaseChange());
        thermo_.readProperty("rho", fluidThermo::LIQUID, thermo_.speciesPhaseChange());

        Sherwood_ =
            SherwoodTypeNames.read(params_.lookup("Sherwood"));

        potential_ =
            potentialTypeNames.read(params_.lookup("potential"));

        diffusion_ =
            diffusionTypeNames.read(params_.lookup("diffusion"));

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
