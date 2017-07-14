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
#include "multiSpeciesFriedlander.H"

#include "makeCondensationEvaporationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{
    makeCondensationEvaporationTypes(multiSpeciesFriedlander, condensationEvaporationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::multiSpeciesFriedlander::multiSpeciesFriedlander
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    condensationEvaporationModel(mesh, aerosol),
    transFunc_(true),
    KelvinEffect_(true)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::condensationEvaporationModels::multiSpeciesFriedlander::~multiSpeciesFriedlander()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesFriedlander::getCondRateListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{
    const dictionary& speciesPhaseChange = thermo_.speciesPhaseChange();
    const label& nSpeciesPhaseChange = thermo_.nSpeciesPhaseChange();

    const dictionary& species = thermo_.species();
    const label& nSpecies = thermo_.nSpecies();

    List<List<scalar> > I(z.size());

    forAll(z, i)
    {
        I[i] = List<scalar>(nSpeciesPhaseChange, 0.0);
    }

    const scalar pi = constant::mathematical::pi;

    const volScalarField& T = thermo_.T();
    const volScalarField& rho = thermo_.rho();
    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const List<scalar>& M = thermo_.M();

    const PtrList<volScalarField>& Y = thermo_.Y();
    const PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<DataEntry<scalar> >& dataEntriesListRho_l =
        thermo_.getProperty("rho", fluidThermo::LIQUID);
    PtrList<DataEntry<scalar> >& dataEntriesListP_s =
        thermo_.getProperty("P_s", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> >& dataEntriesListSigma =
        thermo_.getProperty("sigma", fluidThermo::LIQUID);

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

    // Molecular mass
    List<scalar> m(nSpecies);
    forAll(species, jj)
    {
        m[jj] = 0.001*M[jj]/N_A();
    }

    // Partial molecular volume
    List<scalar> v(nSpeciesPhaseChange);

    // Partial vapor pressure
    List<scalar> p_v(nSpeciesPhaseChange);

    // Saturation vapor pressure
    List<scalar> p_s(nSpeciesPhaseChange);

    // Liquid density of i-species
    List<scalar> rho_l(nSpeciesPhaseChange);

    // Saturation
    List<scalar> S(nSpeciesPhaseChange);

    //Surface tension
    List<scalar> sigma(nSpeciesPhaseChange);

    forAll(speciesPhaseChange, jj)
    {
        // Liquid density
        rho_l[jj] = dataEntriesListRho_l[jj].value(T[jCell]);

        // Partial vapor pressure
        p_v[jj] = max(Y[jj][jCell],0.0)*rho[jCell]*k()*T[jCell]/m[jj];

        // Saturation vapor pressure
        p_s[jj] = dataEntriesListP_s[jj].value(T[jCell]);

        // Saturation
        S[jj] = p_v[jj]/p_s[jj];

        // Partial molecular volume
        v[jj] = m[jj]/rho_l[jj];

        // Surface tension
        sigma[jj] = dataEntriesListSigma[jj].value(T[jCell]);
    }

    scalar sumZm = 0.0;
    scalar sumZ = 0.0;
    scalar sumZrho = 0.0;

    forAll(speciesPhaseChange, jj)
    {
        sumZm += max(Z[jj][jCell],0.0)/m[jj];
        sumZ +=  max(Z[jj][jCell],0.0);
        sumZrho +=  max(Z[jj][jCell],0.0)/stabilise(rho_l[jj], SMALL);
    }

    scalar rho_d = sumZ/stabilise(sumZrho, SMALL);

    scalar sumY = 0.0;
    scalar sumYm = 0.0;

    forAll(species, jj)
    {
        sumY +=  max(Y[jj][jCell],0.0);
        sumYm +=  max(Y[jj][jCell],0.0)/m[jj];
    }

    scalar sigmaDrop = 0.0;
    List<scalar> Wi(nSpeciesPhaseChange);

    forAll(speciesPhaseChange, jj)
    {
        // Mole fraction in the droplet phase
        if (sumZm > 0.0)
        {
            Wi[jj] = (max(Z[jj][jCell],0.0) / m[jj] ) / sumZm;
        }
        else
        {
            Wi[jj] = 0.0;
        }

        // Surface tension of composite droplet
        sigmaDrop += Wi[jj] * sigma[jj];
    }

    scalar mg = sumY/sumYm;

    scalar lambda = sqrt(8.0*k()*T[jCell]/pi/mg)*(4.0*muEff[jCell]/5.0/(p1[jCell]+p0.value()));

    forAll(speciesPhaseChange, j)
    {
        if (Wi[j] > 0.0 )
        {
            forAll(z, i)
            {

                scalar f = 1.0;
                scalar E = 1.0;

                scalar d = pow(z[i]/rho_d * 6.0/pi, 1.0/3.0);

                if (transFunc_ && d > SMALL)
                {
                    f = (1.0 + 2.0*lambda/d)/(1.0 + 5.33*sqr(lambda/d) + 3.42*lambda/d);
                }

                if (KelvinEffect_ && d > SMALL)
                {
                    E = exp(4.0*sigmaDrop*v[j]/(k()*T[jCell]*d));
                }

                // Diffusivity with respect to last species, for now

                scalar D12 = thermo_.getDiffusivity(j, thermo_.nSpecies()-1, T[jCell], p1[jCell]+p0.value());

                scalar Xis = Wi[j]*p_s[j]/(p1[jCell]+p0.value());
                scalar Yis = Xis*m[j]/(Xis*m[j]+(1.0-Xis)*mg);

                I[i][j] = 2.0*pi * d * rho[jCell] * f * D12*Yis*(S[j]-E);
            }
        }
    }

    return I;
}

Foam::List<Foam::List<Foam::scalar> >
Foam::condensationEvaporationModels::multiSpeciesFriedlander::getEtaGammaListCell
(
    const Foam::List<Foam::scalar>& z,
    const Foam::label jCell
)
{
    const dictionary& speciesPhaseChange = thermo_.speciesPhaseChange();
    const label& nSpeciesPhaseChange = thermo_.nSpeciesPhaseChange();

    const dictionary& species = thermo_.species();
    const label& nSpecies = thermo_.nSpecies();

    List<List<scalar> > I(z.size());

    forAll(z, i)
    {
        I[i] = List<scalar>(nSpeciesPhaseChange+1, 0.0);
    }

    const scalar pi = constant::mathematical::pi;

    const volScalarField& T = thermo_.T();
    const volScalarField& rho = thermo_.rho();
    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const List<scalar>& M = thermo_.M();

    const PtrList<volScalarField>& Y = thermo_.Y();
    const PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<DataEntry<scalar> >& dataEntriesListRho_l =
        thermo_.getProperty("rho", fluidThermo::LIQUID);
    PtrList<DataEntry<scalar> >& dataEntriesListP_s =
        thermo_.getProperty("P_s", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> >& dataEntriesListSigma =
        thermo_.getProperty("sigma", fluidThermo::LIQUID);

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

    // Molecular mass
    List<scalar> m(nSpecies);
    forAll(species, jj)
    {
        m[jj] = 0.001*M[jj]/N_A();
    }

    // Partial molecular volume
    List<scalar> v(nSpeciesPhaseChange);

    // Partial vapor pressure
    List<scalar> p_v(nSpeciesPhaseChange);

    // Saturation vapor pressure
    List<scalar> p_s(nSpeciesPhaseChange);

    // Liquid density of i-species
    List<scalar> rho_l(nSpeciesPhaseChange);

    // Saturation
    List<scalar> S(nSpeciesPhaseChange);

    //Surface tension
    List<scalar> sigma(nSpeciesPhaseChange);

    forAll(speciesPhaseChange, jj)
    {
        // Liquid density
        rho_l[jj] = dataEntriesListRho_l[jj].value(T[jCell]);

        // Partial vapor pressure
        p_v[jj] = max(Y[jj][jCell],0.0)*rho[jCell]*k()*T[jCell]/m[jj];

        // Saturation vapor pressure
        p_s[jj] = dataEntriesListP_s[jj].value(T[jCell]);

        // Saturation
        S[jj] = p_v[jj]/p_s[jj];

        // Partial molecular volume
        v[jj] = m[jj]/rho_l[jj];

        // Surface tension
        sigma[jj] = dataEntriesListSigma[jj].value(T[jCell]);
    }

    scalar sumZm = 0.0;
    scalar sumZ = 0.0;
    scalar sumZrho = 0.0;

    forAll(speciesPhaseChange, jj)
    {
        sumZm +=  max(Z[jj][jCell],0.0)/m[jj];
        sumZ +=  max(Z[jj][jCell],0.0);
        sumZrho += max(Z[jj][jCell],0.0)/stabilise(rho_l[jj], SMALL);
    }

    scalar rho_d = sumZ/stabilise(sumZrho, SMALL);

    scalar sumY = 0.0;
    scalar sumYm = 0.0;

    forAll(species, jj)
    {
        sumY += max(Y[jj][jCell],0.0);
        sumYm += max(Y[jj][jCell],0.0)/m[jj];
    }

    scalar sigmaDrop = 0.0;
    List<scalar> Wi(nSpeciesPhaseChange);

    forAll(speciesPhaseChange, jj)
    {
        // Mole fraction in the droplet phase
        if (sumZm > 0.0)
        {
            Wi[jj] = (max(Z[jj][jCell],0.0) / m[jj] ) / sumZm;
        }
        else
        {
            Wi[jj] = 0.0;
        }

        // Surface tension of composite droplet
        sigmaDrop += Wi[jj] * sigma[jj];
    }

    scalar mg = sumY/sumYm;

    scalar lambda = sqrt(8.0*k()*T[jCell]/pi/mg)*(4.0*muEff[jCell]/5.0/(p1[jCell]+p0.value()));

    scalar a = 2.0/3.0/(pow(aerosol_.x()[aerosol_.P()-1],2.0/3.0)-pow(aerosol_.x()[0],2.0/3.0));

    scalar b = - pow(aerosol_.x()[0],2.0/3.0) * a * 3.0/2.0;

    if (sumZ > 0.0)
    {
        forAll(z, i)
        {
            I[i][nSpeciesPhaseChange] = 3.0/2.0 * a * pow(z[i], 2.0/3.0) + b;

            scalar f = 1.0;

            scalar d = pow(z[i]/rho_d * 6.0/pi, 1.0/3.0);

            if (transFunc_ && d > SMALL)
            {
                f = (1.0 + 2.0*lambda/d)/(1.0 + 5.33*sqr(lambda/d) + 3.42*lambda/d);
            }

            forAll(speciesPhaseChange, j)
            {
                scalar E = 1.0;

                if (KelvinEffect_ && d > SMALL)
                {
                    E = exp(4.0*sigmaDrop*v[j]/(k()*T[jCell]*d));
                }

                // Diffusivity with respect to last species, for now

                scalar D12 = thermo_.getDiffusivity(j, thermo_.nSpecies()-1, T[jCell], (p1[jCell]+p0.value()));

                scalar Xis = Wi[j]*p_s[j]/(p1[jCell]+p0.value());
                scalar Yis = Xis*m[j]/(Xis*m[j]+(1.0-Xis)*mg);

                I[i][j] = a * 2.0*pi * pow(6.0/pi/rho_d, 1.0/3.0) * rho[jCell] * f * D12*Yis*(S[j]-E);
            }
        }
    }

    return I;
}

Foam::List<Foam::scalar>
Foam::condensationEvaporationModels::multiSpeciesFriedlander::psiInv
(
    const Foam::List<Foam::scalar>& zeta,
    const Foam::label jCell
)
{
    List<scalar> z(zeta.size(), -1.0);

    scalar a = 2.0/3.0/(pow(aerosol_.x()[aerosol_.P()-1],2.0/3.0)-pow(aerosol_.x()[0],2.0/3.0));

    scalar b = - pow(aerosol_.x()[0],2.0/3.0) * a * 3.0/2.0;

    forAll(zeta, i)
    {
        if(zeta[i] >= 0)
        {
            z[i] = pow((zeta[i]-b)/(3.0*a/2.0), 3.0/2.0);
        }
    }

    return z;
}

bool Foam::condensationEvaporationModels::multiSpeciesFriedlander::read()
{
    if (condensationEvaporationModel::read())
    {
        params_.lookup("transitionFunction") >> transFunc_;
        params_.lookup("KelvinEffect") >> KelvinEffect_;

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
