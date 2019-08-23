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

#include "coupledNucleation.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void coupledNucleation::computeComposition
(
    scalarList& w,
    scalarList& pm,
    scalar& alpha,
    const scalar& f,
    const scalar& Ke,
    const scalar& p,
    const scalarList& pVap,
    const scalarList& pSat,
    const scalarList& D,
    const scalarList& v
) const
{
    #define F (sum(S*Foam::exp(-gamma*beta))-1.0)
    #define dFdb sum(-gamma*S*Foam::exp(-gamma*beta))

    scalar betaOld(0.0);

    const scalar vm(sum(v)/v.size());
    const scalarList gamma(v/vm);

    for (label outerIter = 0; outerIter <= (maxOuterIter_-1); outerIter++)
    {
        // Compute xi

        const scalarList pSurf(f*Ke*w*pSat);

        const scalarList xi(D*Foam::log(max(1.0-sum(pSurf)/p, 0.1)));

        pm = Foam::exp(xi)*pVap;

        const scalarList S(pm/pSat);

        scalar beta(betaOld);

        if (outerIter == 0)
        {
            const scalar maxSat(findMax<scalarList>(S));
            beta = Foam::log(S[maxSat])/gamma[maxSat];
        }

        for (label innerIter = 0; innerIter <= (maxInnerIter_-1); innerIter++)
        {
            const scalar betaNew = beta - F/dFdb;

            if (mag(betaNew-beta) < TOL_)
            {
                beta = betaNew;

                break;
            }
            else if (innerIter == (maxInnerIter_-1))
            {
                WarningInFunction
                    << "The model did not converge in the inner loop. "
                    << " tolerance = " << mag(betaNew-beta)
                    << endl;

                beta = betaNew;

                break;
            }

            beta = betaNew;
        }

        w = S*Foam::exp(-beta*gamma);

        if (mag(betaOld-beta) < TOL_)
        {
            alpha = beta/vm;

            break;
        }
        else if (outerIter == (maxOuterIter_-1))
        {
            WarningInFunction
                << "The model did not converge in the outer loop."
                << " tolerance = " << mag(betaOld-beta)
                << endl;

            alpha = beta/vm;

            break;
        }

        betaOld = beta;
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledNucleation, 0);
addToRunTimeSelectionTable(nucleationModel, coupledNucleation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledNucleation::coupledNucleation
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    nucleationModel(type(), aerosol, dict),
    TOL_(readScalar(dict.lookup("tolerance"))),
    maxInnerIter_(dict.lookupOrDefault<scalar>("maxInnerIter", 100)),
    maxOuterIter_(dict.lookupOrDefault<scalar>("maxOuterIter", 100))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coupledNucleation::~coupledNucleation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

nucData coupledNucleation::rate
(
    const scalar& p,
    const scalar& T,
    const scalarList& Y,
    const scalarList& pSat,
    const scalarList& D,
    const scalarList& rhoDisp,
    const scalarList& sigma
) const
{
    const scalar pi = constant::mathematical::pi;
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar kB = constant::physicoChemical::k.value();

    aerosolThermo& thermo = aerosol_.thermo();

    const speciesTable& activeSpecies = thermo.activeSpecies();

    const basicSpecieMixture& compCont = thermo.thermoCont().composition();

    // Prepare data

    const scalarList Ya(Y, thermo.activeSpeciesMap());
    const scalarList Yia(Y, thermo.inactiveSpeciesMap());

    scalarList W(Y.size(), 0.0);

    forAll(Y, j)
    {
        W[j] = compCont.W(j);
    }

    const scalarList Wa(W, thermo.activeSpeciesMap());

    const scalar sumY(min(sum(Y), 1.0));
    const scalar sumYa(min(sum(Ya), 1.0));
    const scalar sumYia(min(sum(Yia), 1.0));

    const scalarList y(Y/sumY);
    const scalarList x(y/W/sum(y/W));

    const scalarList xa(x, thermo.activeSpeciesMap());

    const scalarList pVap(p*xa);

    scalarList S(pVap/pSat);
    const label maxSat(findMax<scalarList>(S));

    // Create output data

    nucData data(activeSpecies.size());

    // Check if we have an adequate mixture

    if
    (
        sumYa > SMALL
     && sumYia > SMALL
     && maxSat > -1
     && S[maxSat] > (1.0+SMALL)
    )
    {
        data.active() = true;

        const scalarList m(0.001*Wa/NA);
        const scalarList v(m/rhoDisp);

        // Compute diffusivities

        const scalarList xia(x, thermo.inactiveSpeciesMap());
        const scalarList Dia(D, thermo.inactiveSpeciesMap());
        const scalarList Da(D, thermo.activeSpeciesMap());

        const scalar DiaMean(sum(Dia*xia)/sum(xia));

        // Kelvin and activity coefficients to unity, for now

        scalar Ke = 1.0;
        scalar f = 1.0;

        // Compute the critical cluster composition

        scalarList w(activeSpecies.size(), 1.0/scalar(activeSpecies.size()));
        scalarList pm(activeSpecies.size(), 0.0);
        scalar alpha(0.0);

        computeComposition
        (
            w,
            pm,
            alpha,
            f,
            Ke,
            p,
            pVap,
            pSat,
            DiaMean/Da,
            v
        );

        S = pm/pSat;

        // Critical cluster properties

        const scalar kBT(kB*T);
        const scalar vmc(sum(w*v));
        const scalar sigmac(sum(w*sigma));
        const scalar rc(2.0*sigmac/(alpha*kBT));
        const scalar vc(pi/6.0*Foam::pow(rc*2.0,3.0));
        const scalar Nc(vc/vmc);

        if (Nc > 1.0)
        {
            const scalar Gc(4.0/3.0*pi*Foam::sqr(rc)*sigmac);
            const scalar mc(Nc*sum(w*m));

            const scalarList s
            (
                Foam::pow(v,2.0/3.0)
              * Foam::pow(36.0*pi,1.0/3.0)
            );

            // Equilibrium critical cluster distribution

            scalar zc(Foam::exp(-Gc/kBT));

            forAll(activeSpecies, j)
            {
                zc *= Foam::pow(pSat[j]/kBT*Foam::exp(s[j]*sigma[j]/kBT),w[j]);
            }

            // Zeldovich factor

            const scalar Zec
            (
                Foam::pow
                (
                    sigmac
                  * Foam::sqr(vmc)
                  / (4.0*kBT*Foam::sqr(pi)*Foam::pow(rc,4.0)),
                    1.0-scalar(activeSpecies.size())/2.0
                )
            );

            // Cluster composition

            const scalarList z(w*Wa/sum(w*Wa));

            // Growth rate

            const scalarList K
            (
                pVap/kBT*Foam::pow(3.0/(4.0*pi),1.0/6.0)*Foam::sqrt(6.0*kBT)
              * Foam::sqrt(1.0/m+1.0/mc)
              * Foam::sqr(Foam::pow(v,1.0/3.0)+Foam::pow(vc,1.0/3.0))
            );

            const scalar Kc
            (
                sum(Foam::sqr(w))
              / sum(Foam::sqr(w)/max(K,VSMALL))
            );

            // Set the nucleation data

            data.J() = Kc*Zec*zc;
            data.s() = mc;
            data.z() = z;
        }
    }

    return data;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
