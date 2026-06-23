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

#include "coupledCondensation.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledCondensation, 0);
addToRunTimeSelectionTable(condensationModel, coupledCondensation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledCondensation::coupledCondensation
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    condensationModel(type(), aerosol, dict),
    KelvinEffect_(dict.lookupOrDefault<Switch>("KelvinEffect", false)),
    FuchsCorrection_(dict.lookupOrDefault<Switch>("FuchsCorrection", false)),
    DropletTemperatureCorrection_(dict.lookupOrDefault<Switch>("DropletTemperatureCorrection",false)), 

    
    SR_(dict.lookupOrDefault<scalar>("SR", 0.0)),   //Default is 0
    soluteName_(dict.lookupOrDefault<word>("solute", "none")),
    soluteLimit_(dict.lookupOrDefault<scalar>("soluteLimit", -1.0))
{
    if (DropletTemperatureCorrection_ && !dict.found("SR"))
    {
        FatalErrorInFunction
            << "DropletTemperatureCorrection is enabled but "
            << "SR is not specified."
            << nl
            << "Please provide SR in the dictionary."
            << exit(FatalError);
    }

    if (DropletTemperatureCorrection_ && SR_ < 0.0)
    {
        FatalErrorInFunction
            << "SR must be >= 0.0. Current value: "
            << SR_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coupledCondensation::~coupledCondensation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

conData coupledCondensation::rate
(
    const scalar d,
    const scalar p,
    const scalar T,
    const scalarList& Y,
    const scalarList& Z,
    const scalarList& pSat,
    const scalarList& gamma,
    const scalarList& D,
    const scalarList& rhoCont,
    const scalarList& rhoDisp,
    const scalarList& sigma
) const
{
    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();
    const scalar NA = constant::physicoChemical::NA.value();

    aerosolThermo& thermo = aerosol_.thermo();
    
    // Search if there is solute
    static label soluteDispIndex = -2;  // -2 not searched, -1 not found
    
        if (soluteDispIndex == -2)
    {
        soluteDispIndex = -1;

        if (soluteName_ != "none")
        {
            const speciesTable& dispSpec =
                thermo.thermoDisp().composition().species();

            forAll(dispSpec, i)
            {
                if (dispSpec[i] == soluteName_)
                {
                    soluteDispIndex = i;
                    break;
                }
            }

            if (soluteDispIndex >= 0)
            {
                Info<< "coupledCondensation: Solute '" << soluteName_
                    << "' found in dispersed phase at index "
                    << soluteDispIndex << endl;
            }
            else
            {
                Warning<< "coupledCondensation: Solute '" << soluteName_
                       << "' NOT found in dispersed phase. Solute effects disabled."
                       << endl;
            }
        }
    }


    rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

    const basicSpecieMixture& compCont = thermoCont.composition();

    const speciesTable& activeSpecies = thermo.activeSpecies();

    const scalarList Ya(Y, thermo.activeSpeciesMap());

    const scalar sumY(min(sum(Y), 1.0));
    const scalar sumYa(min(sum(Ya), 1.0));
    const scalar sumZ(min(sum(Z), 1.0));

    conData data(activeSpecies.size());
    
    // -------------------------------------------------
    // Solute limiting condition (if enabled)
    // -------------------------------------------------
    if (soluteDispIndex >= 0 && soluteLimit_ > 0.0 && sumZ > 1E-30 && (sumY - sumYa) > 0.0)
    {
        const scalar Zsolute = Z[soluteDispIndex] / sumZ;

        if (Zsolute >= soluteLimit_)
        {
            data.active() = false;

            forAll(data.source(), j)
            {
                data.source()[j] = 0.0;
                data.sink()[j]   = 0.0;
            }

            return data;
        }
    }

    // -------------------------------------------------
    // Normal condensation

    // Check if we have an adequate mixture

    if (sumZ > 1E-30 && (sumY-sumYa) > 0.0)
    
    {   

        data.active() = true;

        scalarList W(Y.size(), 0.0);

        forAll(Y, j)
        {
            W[j] = compCont.W(j);
        }

        const scalarList Wd(W, thermo.dispSpeciesMap());

        // Compute fractions w.r.t. the dispersed phase

        const scalarList z(Z/sumZ);
        const scalarList w(z/Wd/sum(z/Wd));

        // Compute fractions w.r.t. the continuous phase

        const scalarList y(Y/sumY);
        const scalarList x(y/W/sum(y/W));


        // Effect of droplet temperature

	scalar dTdrop = 0.0;

	if (DropletTemperatureCorrection_)
	{
    		dTdrop = (
          ((6.65 + 0.345*(T-273.15)
           + 0.0031*sqr(T-273.15))* (SR_ - 1.0))/(1.0 + (0.082 + 0.00782*(T-273.15))*SR_)
    	                 );
	}

        const scalar DropletTemp = T + dTdrop;
         
        // Kelvin effect factor

        scalar Ke = 1.0;

        if (KelvinEffect_)
        {
            const scalarList Md(0.001*Wd/NA);
            const scalar rhol = 1.0/sum(z/rhoDisp);
            const scalar vl = sum(w*Md)/rhol;
            const scalar sigmal = sum(w*sigma);

            Ke = exp(4.0*sigmal*vl/(kB*DropletTemp*d));
        }

        // Set Fuchs & Sutugin correction factor

        scalar beta = 1.0;

        if (FuchsCorrection_)
        {	
       	// Compute lambda
       
       	const label jInert(thermoCont.species()[aerosol_.thermo().inertSpecie()]);
       	const scalar mu(compCont.mu(jInert, p, T));
       	const scalar mg(0.001*W[jInert]/NA);
       
       	const scalar lambda
    			(
        			Foam::sqrt(8.0*kB*T/(pi*mg)) * 4.0/5.0*mu/p
    			);
            
       	const scalar Kn = 2.0 * lambda / max(d, aerosol_.dMin());
       	
       	beta = (1.0+ Kn)/(1.0 +1.71 * Kn + 1.333 * Kn * Kn);
       
        }

        
        // Compute pressures

        
        const scalarList pVap(p*x);
        const scalarList pVapOverY(p/W/sum(Y/W));
        
        scalarList pSurf(activeSpecies.size(), 0.0);
        scalarList pSurfOverZ(activeSpecies.size(), 0.0);
        
        
        if (DropletTemperatureCorrection_)
	{
    		// Compute saturation vapour pressure at droplet temperature // Pa (eq. 13.2 Hinds, 1999)
    		const scalar pSat_DropletTemp =
        	1e3*Foam::exp(16.7 - (4060.0/(DropletTemp - 37.0)));

          pSurf = gamma*Ke*pSat_DropletTemp*w;

          pSurfOverZ = gamma*Ke*pSat_DropletTemp/Wd/sum(Z/Wd);
        }
        else
        {
          pSurf = gamma*Ke*pSat*w;

          pSurfOverZ = gamma*Ke*pSat/Wd/sum(Z/Wd);
        }
        

        // Compute xi

        // Compute diffusivities

        const scalarList xia(x, thermo.inactiveSpeciesMap());
        const scalarList Dia(D, thermo.inactiveSpeciesMap());
        const scalarList Da(D, thermo.activeSpeciesMap());

        const scalar DiaMean(sum(Dia*xia)/sum(xia));

        const scalarList pSurfa(pSurf, thermo.activeSpeciesMap());
        const scalarList pVapa(pVap, thermo.activeSpeciesMap());

        const scalarList xi
        (
            DiaMean/Da
          * Foam::log(max((p-sum(pSurfa))/(p-sum(pVapa)), 0.1))
        );

        // Compute condensation rates

        forAll(activeSpecies, j)
        {
            scalar func(0.0);

            if (mag(xi[j]) < 1E-10)
            {
                func = -1.0 + 0.5*xi[j];
            }
            else
            {
                func = xi[j] / (1.0-Foam::exp(xi[j]));
            }

            const scalar c(-2.0*pi*beta*D[j]*rhoCont[j]*func/p);

            data.source()[j] = c*Foam::exp(xi[j])*pVapOverY[j];
            data.sink()[j] = c*pSurfOverZ[j];
        }
    }

    return data;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
