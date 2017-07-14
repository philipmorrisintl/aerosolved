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
#include "multiSpeciesCNT.H"

#include "makeNucleationTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{
    makeNucleationTypes(multiSpeciesCNT, nucleationModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationModels::multiSpeciesCNT::multiSpeciesCNT
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

Foam::nucleationModels::multiSpeciesCNT::~multiSpeciesCNT()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField>&
Foam::nucleationModels::multiSpeciesCNT::getNucFields()
{
    const scalar pi = constant::mathematical::pi;

    const volScalarField& T = thermo_.T();
    const volScalarField& rho = thermo_.rho();

    const dictionary speciesPhaseChange = thermo_.speciesPhaseChange();
    const label n = thermo_.nSpeciesPhaseChange();

    const dictionary& species = thermo_.species();
    const label& nSpecies = thermo_.nSpecies();

    const List<scalar>& M = thermo_.M();
    const List<scalar>& Tc = thermo_.Tc();

    const PtrList<volScalarField>& Y = thermo_.Y();

    PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo_.getProperty("rho", fluidThermo::LIQUID);
    PtrList<DataEntry<scalar> > dataEntriesListP_s =
        thermo_.getProperty("P_s", fluidThermo::VAPOR);
    PtrList<DataEntry<scalar> > dataEntriesListSigma =
        thermo_.getProperty("sigma", fluidThermo::LIQUID);

    // Molecular mass
    List<scalar> m(nSpecies);
    forAll(species, i)
    {
       m[i] = 0.001*M[i]/N_A();
    }

    // Partial molecular volume
    List<scalar> v(n);
    List<scalar> fi(n);

    // Partial vapor pressure
    List<scalar> p_v(n);

    // Saturation vapor pressure
    List<scalar> p_s(n);

    // Saturation
    List<scalar> S(n);

    // Liquid density of i-species
    List<scalar> rho_l(n);

    // Number of molecules of i-species in the critical cluster
    List<scalar> Ni(n, 0.0);

    //Surface tension
    List<scalar> sigma(n);

    forAll(mesh_.C(), jCell)
    {
        scalar maxS = 0.0;
        scalar sumYPhaseChange = 0.0;

        // 0-1 list of if T is below the Tc of the species
        List<scalar> belowTc(n);

        forAll(speciesPhaseChange, i)
        {
            // Check if T is below the Tc of the species
            belowTc[i] = (T[jCell] <= Tc[i]) ? 1.0 : 0.0;

            // Liquid density
            rho_l[i] = dataEntriesListRho_l[i].value(min(T[jCell],Tc[i]));

            // Partial molecular volume
            v[i] = m[i]/rho_l[i];

            // Partial vapor pressure
            p_v[i] = max(Y[i][jCell],0.0)*rho[jCell]*k()*T[jCell]/m[i];

            // Saturation vapor pressure
            p_s[i] = dataEntriesListP_s[i].value(T[jCell]);

            // Saturation
            S[i] = p_v[i]/p_s[i];
            fi[i] = v[i]/v[0];

            // Max saturation
            maxS = max(maxS,S[i]);

            // Surface tension
            sigma[i] = dataEntriesListSigma[i].value(min(T[jCell],Tc[i]));

            sumYPhaseChange += belowTc[i]*max(Y[i][jCell],0.0);
        }

        // Compute minimal diameter based on the min(v) - min molecular volume
        scalar minDiameter = 2.0*pow(3.0*min(v)/4.0/pi,1./3.);

        SJDnuc_[n][jCell] = 0.0;
        SJDnuc_[n+1][jCell] = 0.0;

        if ( sumYPhaseChange > TOL() && maxS > (1.0+TOL()) )
        {
            scalar ksi = -log(1.0/maxS);
            scalar ksistart = 0.0;

            // Newton method for alpha
            int iter = 0;
            while ( abs((ksi - ksistart)/ksi) > TOL() )
            {
                iter++;
                ksistart = ksi;
                scalar fSum  = 0.0;
                scalar dfSum = 0.0;
                forAll(speciesPhaseChange, i)
                {
                    fSum  += S[i]*exp( - ksi * fi[i]  );
                    dfSum -=  fi[i]*S[i]*exp( - ksi * fi[i] );
                }
                ksi = ksistart - (fSum - 1.0 ) / dfSum;
            }
            scalar alpha = ksi/v[0];

            // Mole fraction in critical cluster
            List<scalar> w(n);

            scalar wSum = 0.0;
            forAll(speciesPhaseChange, i)
            {
                w[i] = S[i]*exp( - alpha * v[i]);
                wSum += w[i];
            }

            scalar vSum = 0.0;
            scalar sigmaSum = 0.0;
            forAll(speciesPhaseChange, i)
            {
                w[i] /= wSum;
                vSum += w[i]*v[i];
                sigmaSum += w[i]*sigma[i];
            }

            // Diameter of critical cluster
            scalar diameter = 4.0*sigmaSum/(k()*T[jCell]*alpha);

            // Gibbs free energy
            scalar deltaG = pi*sqr(diameter)*sigmaSum/3.0;

            // Total number of molecules
            scalar Ntot = pi*pow(diameter,3.0)/vSum/6.0;

            // Total mass of the cluster
            scalar mSum = 0.0;
            scalar nnuc = 0.0;

            // Heaviside function result H(w[i])
            List<scalar> H(n);

            //Condensation rate of i-species
            List<scalar> Kii(n);

            forAll(speciesPhaseChange, i)
            {
                Ni[i] = Ntot * w[i];

                H[i] = (w[i] > 0.0) ? 1.0 : 0.0;

                nnuc += H[i];
                mSum += Ni[i] * m[i];
            }

            // Equilibrium concentration of critical cluster
            scalar ceq = 0.0;

            scalar sqrNSum = 0.0;
            scalar sqrNSumKii = 0.0;
            scalar fac = pow(3.0/pi/4.0,1.0/6.0)*sqrt(6.0*k()*T[jCell])/k()/T[jCell];

            scalar PiTerm = 1.0;
            forAll(speciesPhaseChange, i)
            {
                scalar simon = pow(36.0*pi,1.0/3.0)*pow(v[i],2.0/3.0);
                PiTerm *= pow(p_s[i]/k()/T[jCell]*exp(simon*sigma[i]/k()/T[jCell]), w[i]);
            }

            forAll(speciesPhaseChange, i)
            {
                // Wilemski correction
                ceq = exp( - deltaG/k()/T[jCell] ) * PiTerm;

                Kii[i] = p_v[i] * fac * sqrt(1./m[i] + 1.0/mSum) *
                sqr(pow(m[i]/rho_l[i],1.0/3.0) + 0.5*diameter*pow(4.0*pi/3.0,1.0/3.0));

                sqrNSum += Ni[i]*Ni[i];
                sqrNSumKii += (Kii[i] != 0.0) ? Ni[i]*Ni[i]/Kii[i] : 0.0;
            }

            // Zeldovich factor
            scalar Ze = 0.0;
            if (diameter > minDiameter )
            {
                Ze = pow(sigmaSum*vSum*vSum/(k()*T[jCell]*sqr(pi)/4.0*pow(diameter,4.0)),(1.0-nnuc/2.0));
            }

            // Average growth rate
            scalar Rav = (sqrNSumKii > 0.0) ? sqrNSum/sqrNSumKii : 0.0;

            // Nucleation rate
            SJDnuc_[n][jCell] = Rav * Ze * ceq;

            // Critical size
            SJDnuc_[n+1][jCell] = mSum*2.0;
        }

        // Condensation by nucleation
        forAll(speciesPhaseChange, j)
        {
            SJDnuc_[j][jCell] = SJDnuc_[n][jCell] * 2.0 * Ni[j] * m[j];
        }
    }

    return SJDnuc_;
}

bool Foam::nucleationModels::multiSpeciesCNT::read()
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

inline const Foam::scalar& Foam::nucleationModels::multiSpeciesCNT::TOL() const
{
    return TOL_;
}

// ************************************************************************* //
