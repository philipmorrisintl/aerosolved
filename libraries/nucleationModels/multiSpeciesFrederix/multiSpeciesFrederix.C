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

#include "nucleationModel.H"
#include "multiSpeciesFrederix.H"

#include "makeNucleationTypes.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{
    makeNucleationTypes(multiSpeciesFrederix, nucleationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationModels::multiSpeciesFrederix::multiSpeciesFrederix
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    nucleationModel(mesh, aerosol)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationModels::multiSpeciesFrederix::~multiSpeciesFrederix()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField>&
Foam::nucleationModels::multiSpeciesFrederix::getNucFields()
{
    const scalar pi = constant::mathematical::pi;
    const scalar kB = k();
    const scalar NA = N_A();

    const label N = thermo_.nSpecies();
    const label Nps = thermo_.nSpeciesPhaseChange();

    const dictionary& setS = thermo_.species();
    const dictionary& setSps = thermo_.speciesPhaseChange();

    const volScalarField& T = thermo_.T();
    const volScalarField& p1 = thermo_.p1();
    const scalar p0 = thermo_.p0().value();

    const List<scalar>& M = thermo_.M();
    const List<scalar>& Tc = thermo_.Tc();

    const PtrList<volScalarField>& Y = thermo_.Y();

    PtrList<DataEntry<scalar> > dataEntriesListP_s =
        thermo_.getProperty("P_s", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo_.getProperty("rho", fluidThermo::LIQUID);
    PtrList<DataEntry<scalar> > dataEntriesListSigma =
        thermo_.getProperty("sigma", fluidThermo::LIQUID);

    // Molecular mass

    scalarList m(N, 0.0);

    forAll(setS, j)
    {
        m[j] = 1E-3*M[j]/NA;
    }

    forAll(mesh_.C(), jCell)
    {
        const scalar kT = kB*T[jCell];

        // Vapor mixture mean molar weight

        scalar sumYj(0.0);
        scalar sumYjOverMj(0.0);

        forAll(setS, j)
        {
            sumYj += max(Y[j][jCell],0.0);
            sumYjOverMj += max(Y[j][jCell],0.0)/M[j];
        }

        const scalar Mv = sumYj/sumYjOverMj;

        // Prepate fields

        scalarList S(Nps, 0.0);
        scalarList pv(Nps, 0.0);
        scalarList ps(Nps, 0.0);
        scalarList rhol(Nps, 0.0);
        scalarList v(Nps, 0.0);
        scalarList sigma(Nps, 0.0);

        scalar sumYps(0.0);

        // 0-1 list of if T is below the Tc of the species
        List<scalar> belowTc(Nps);

        forAll(setSps, j)
        {
            // Check if T is below the Tc of the species
            belowTc[j] = (T[jCell] <= Tc[j]) ? 1 : 0;

            pv[j] = max(Y[j][jCell],0.0) * Mv / M[j]
                  * (p0 + p1[jCell]);

            ps[j] = dataEntriesListP_s[j].value(T[jCell]);

            rhol[j] = dataEntriesListRho_l[j].value(min(T[jCell],Tc[j]));

            sigma[j] = dataEntriesListSigma[j].value(min(T[jCell],Tc[j]));

            S[j] = pv[j]/ps[j];

            v[j] = m[j]/rhol[j];

            sumYps += belowTc[j]*max(Y[j][jCell],0.0);
        }

        const label l = findMax<scalarList>(S);

        SJDnuc_[Nps][jCell] = 0.0;
        SJDnuc_[Nps+1][jCell] = 0.0;

        if (sumYps > TOL_ && l > -1 && S[l] > (1.0+TOL_))
        {
            // Solve for beta

            const scalar vm = sum(v)/v.size();

            const scalarList gamma = v/vm;

            #define F (sum(S*exp(-gamma*b))-1.0)
            #define dFdb sum(-gamma*S*exp(-gamma*b))

            scalar b = log(S[l])/gamma[l];

            for (label i = 0; i <= 9; i++)
            {
                const scalar bNew = b - F/dFdb;

                if (mag(bNew-b) < TOL_)
                {
                    b = bNew;

                    break;
                }
                else if (i == 9)
                {
                    FatalErrorIn
                    (
                        "Foam::nucleationModels::multiSpeciesFrederix::getNucFields()"
                    )   << "The model did not converge in cell " << jCell << nl
                        << exit(FatalError);
                }
                else
                {
                    b = bNew;
                }
            }

            // Critical cluster properties

            const scalarList omegac = S*exp(-b*gamma);

            const scalar vmc = sum(omegac*v);

            const scalar sigmac = sum(omegac*sigma);

            const scalar rc = 2.0*sigmac / (kT * b/vm);

            const scalar vc = pi/6.0 * pow(rc*2.0, 3.0);

            const scalar Gc = 4.0/3.0 * pi * sqr(rc) * sigmac;

            const scalarList s = pow(v, 2.0/3.0) * pow(36.0*pi, 1.0/3.0);

            scalar Mc = 0.0;

            forAll(setSps, j)
            {
                Mc += omegac[j]*M[j];
            }

            scalarList wc(Nps, 0.0);

            forAll(setSps, j)
            {
                wc[j] = omegac[j]*M[j]/Mc;
            }

            const scalar Nc = vc/vmc;

            const scalar mc = Nc*sum(omegac*m);

            if (Nc > 1.0)
            {
                // Equilibrium critical cluster distribution

                scalar zc = exp(-Gc/kT);

                forAll(setSps, j)
                {
                    zc *= pow
                        (
                            ps[j]/kT
                          * exp(s[j] * sigma[j] / kT),
                            omegac[j]
                        );
                }

                // Zeldovich factor

                const scalar Zec = pow
                    (
                        sigmac*sqr(vmc) / (kT * 4.0 * sqr(pi) * pow(rc, 4.0)),
                        1.0 - scalar(Nps)/2.0
                    );

                // Growth rate

                const scalarList K = pv/kT * pow(3.0/(4.0*pi), 1.0/6.0) * sqrt(6.0*kT)
                      * sqrt(1.0/m + 1/mc)
                      * sqr(pow(v, 1.0/3.0) + pow(vc, 1.0/3.0));

                const scalar Kc = sum(sqr(omegac)) / sum(sqr(omegac)/stabilise(K,VSMALL));

                // Nucleation rate

                SJDnuc_[Nps][jCell] = Kc * Zec * zc;

                // Mass transfer rates

                forAll(setSps, j)
                {
                    SJDnuc_[j][jCell] = SJDnuc_[Nps][jCell] * mc * wc[j];
                }

                // Critical mass

                SJDnuc_[Nps+1][jCell] = mc;
            }
            else
            {
                SJDnuc_[Nps][jCell] = 0.0;
                SJDnuc_[Nps+1][jCell] = 0.0;

                forAll(setSps, j)
                {
                    SJDnuc_[j][jCell] = 0.0;
                }
            }
        }
    }

    return SJDnuc_;
}

bool Foam::nucleationModels::multiSpeciesFrederix::read()
{
    if (nucleationModel::read())
    {
        // Read nucleation model coefficients

        coeffs_.lookup("TOL") >> TOL_;

        // Read species data

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
