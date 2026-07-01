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

#include "addToRunTimeSelectionTable.H"
#include "fixedSectional.H"
#include "fv.H"
#include "fvOptions.H"
#include "constants.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aerosolModels
{
	inline scalar computeFractalDiameter(scalar ds, scalar N, scalar Df)
	{
    		return ds * pow(N, 1.0 / Df);
	}

	inline scalar computeMobilityDiameter(scalar ds, scalar N, scalar Df, scalar dfractal)
	{
    		scalar t1 = (1.0 / 3.0) * Df * pow(N, (2.0 / Df));
    		scalar t2 = 1.0 + (2.0 / 3.0) * (N - 1);
    		scalar t3 = min(t1, t2);
    		scalar t4 = pow(N, (2.0 / 3.0));
    		scalar t5 = max(t3, t4);
    		scalar t6 = ds * pow(t5, 0.5);
    		scalar t7 = dfractal * pow(0.5 * (Df - 1.0), 0.7);
    		scalar t8 = dfractal / (log(dfractal / ds) + 1.0);
    		return max(max(t6, t7), t8);
	}

} 
} 

namespace Foam
{
namespace aerosolModels
{
    defineTypeNameAndDebug(fixedSectional, 0);
    addToRunTimeSelectionTable(aerosolModel, fixedSectional, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::aerosolModels::fixedSectional::updateDrift()
{
    surfaceScalarField& phiInertial = drift_->phiInertial();
    volScalarField& DDisp = drift_->DDisp();
    surfaceScalarField& phiCorrDiff = drift_->phiCorrDiff();
    volVectorField& Ur = drift_->Ur();

    phiInertial *= 0.0;
    DDisp *= 0.0;
    phiCorrDiff *= 0.0;
    Ur *= 0.0;

    const volScalarField alpha(system_->alpha());
    const volScalarField& rho = this->rho();

    // Reconstruct the inertial drift of mass fraction from the sectional
    // inertial drifts

    if (this->drift().dispInertialDrift().type() != "none")
    {
        volScalarField alpha
        (
            IOobject
            (
                "alpha",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("alpha", dimless, 0)
        );

        forAll(system_->distribution(), i)
        {
            const word sectionName(system_->distribution()[i].sectionName());

            const volScalarField& Mi = system_->distribution()[i].M();

            const volScalarField d(system_->d(i));

            const volVectorField& Vi =
                this->drift().dispInertialDrift().V(d, sectionName);

            surfaceScalarField& phiInertiali =
                system_->distribution()[i].phiInertial();

            phiInertiali = fvc::flux(rho*Vi);

            const volScalarField alphai
            (
                system_->distribution()[i].xd()
              * Mi
            );

            Ur += Vi*alphai;
            alpha += alphai;
        }

        Ur /= max(alpha, residualAlpha_);
        phiInertial = fvc::flux(rho*Ur);
    }
    else
    {
        forAll(system_->distribution(), i)
        {
            system_->distribution()[i].phiInertial() *= 0.0;
        }
    }

    // Reconstruct the diffusivity of mass fraction from the sectional
    // diffusivities

    if (this->drift().dispDiff().type() != "none")
    {
        volScalarField alpha
        (
            IOobject
            (
                "alpha",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("alpha", dimless, 0)
        );

        forAll(system_->distribution(), i)
        {
            const volScalarField& Mi = system_->distribution()[i].M();

            volScalarField& Di = system_->distribution()[i].D();

            Di = this->drift().dispDiff().D(system_->d(i));

            const volScalarField alphai
            (
                system_->distribution()[i].xd()
              * Mi
            );

            DDisp += Di*alphai;
            alpha += alphai;
        }

        DDisp /= max(alpha, residualAlpha_);

        forAll(system_->distribution(), i)
        {
            const volScalarField& Di = system_->distribution()[i].D();
            const volScalarField& Mi = system_->distribution()[i].M();

            phiCorrDiff +=
                linearInterpolate(rho*Di)*system_->distribution()[i].xd()
              * fvc::snGrad(Mi/max(alpha, residualAlpha_))
              * mesh_.magSf();
        }

        phiCorrDiff += linearInterpolate(DDisp)*fvc::snGrad(rho)*mesh_.magSf();
    }
    else
    {
        forAll(system_->distribution(), i)
        {
            system_->distribution()[i].D() *= 0.0;
        }
    }
}

void Foam::aerosolModels::fixedSectional::solveSpatial()
{
    Info<<"fixedSectional: solving spatial step" << endl;

    const volScalarField& rho = this->rho();

    const surfaceScalarField& phi = this->phi();
    const surfaceScalarField& phiCorr = drift_->phiCorr();

    const fv::convectionScheme<scalar>& mvPhi = drift_->mvPhi();
    const fv::convectionScheme<scalar>& mvPhiCorr = drift_->mvPhiCorr();

    // Solve the system of equations

    const volScalarField mut(turbulence().mut());

    forAll(system_->distribution(), i)
    {
        volScalarField& Mi = system_->distribution()[i].M();

        const volScalarField& Di = system_->distribution()[i].D();

        // Correction term to account for the fact that diffusion is
        // proportional to number concentration per unit volume, not per unit
        // mass

        const surfaceScalarField phiCorrDiff
        (
            linearInterpolate(Di)*fvc::snGrad(rho)*mesh_.magSf()
        );

        // TODO: is creation here really a good idea?
        fv::options& fvOptions(fv::options::New(mesh_));

        // Dispersed inertial drift and dispersed diffusive drift

        const surfaceScalarField& phiInertiali =
            system_->distribution()[i].phiInertial();

        fvScalarMatrix drift
        (
            fvm::div(phiInertiali, Mi, "div(mvConv)")
          - fvm::laplacian(mut+rho*Di, Mi, "laplacian(mut+rho*D,M)")
          - fvm::div(phiCorrDiff, Mi, "div(mvConv)")
        );

        // Number concentration equation

        fvScalarMatrix MEqn
        (
            fvm::ddt(rho, Mi)
          + mvPhi.fvmDiv(phi, Mi)
          + drift
          ==
            fvOptions(rho, Mi)
          + mvPhiCorr.fvmDiv(phiCorr, Mi)
        );

        MEqn.relax();

        fvOptions.constrain(MEqn);

        MEqn.solve(mesh_.solver("M"));

        fvOptions.correct(Mi);

        Mi.max(0.0);

        phiEff_[i] = MEqn.flux();
    }

    if (rescale_)
    {
        system_->rescale();
    }

    system_->collect();
}

void Foam::aerosolModels::fixedSectional::solveInternal()
{
    clearRates();

    if
    (
        nucleation_->modelType() == "none"
     && condensation_->modelType() == "none"
     && coalescence_->modelType() == "none"
    )
    {
        return;
    }

    Info<<"fixedSectional: solving internal step" << endl;

    //const scalar pi = constant::mathematical::pi;

    const speciesTable& activeSpecies = thermo_.activeSpecies();
    const speciesTable& contSpecies = thermo_.contSpecies();

    const scalarField& p = thermo_.p().field();
    const scalarField& T = thermo_.T().field();
    const scalarField& rho = this->rho().field();

    //    scalarField rDeltaT(rho.size(),1/mesh_.time().deltaTValue());
    const scalarField &rDeltaT=getRDeltaT();

    PtrList<volScalarField>& Y = thermo_.Y();
    PtrList<volScalarField>& Z = thermo_.Z();

    PtrList<scalarField> pSat(thermo_.pSat(activeSpecies));
    PtrList<scalarField> gamma(thermo_.activity().gamma(activeSpecies));
    PtrList<scalarField> D(thermo_.diffusivity().Deff());
    PtrList<scalarField> sigma(thermo_.sigma(activeSpecies));
    PtrList<scalarField> rhoDisp(thermo_.rhoDisp(activeSpecies));

    const sectionalDistribution& dist = system_->distribution();
    sectionalInterpolation& interp = system_->interpolation();

    PtrList<section>& sections = system_->distribution().sections();

    const scalarField dcm(this->meanDiameter(1,0));
    const scalarField rhol(thermo_.thermoDisp().rho());

    // Nucleation

    if (nucleation_->modelType() != "none")
    {
        forAll(rho, celli)
        {
            const nucData ndata
            (
                nucleation_->rate
                (
                    p[celli],
                    T[celli],
                    entryList(Y,celli),
                    entryList(pSat,celli),
                    entryList(gamma,celli),
                    entryList(D,celli),
                    entryList(rhoDisp,celli),
                    entryList(sigma,celli)
                )
            );

            if (ndata.active())
            {
                J_.field()[celli] = ndata.J();

                interp.addToM
                (
                    ndata.s(),
                    ndata.J()/rho[celli]/rDeltaT[celli],
                    celli
                );

                const scalar Inuc(ndata.s()*ndata.J()/rho[celli]);

                forAll(activeSpecies, j)
                {
                    const scalar dZj
                    (
                        min(Inuc*ndata.z()[j], Y[j][celli])
                      / rDeltaT[celli]
                    );

                    Z[j][celli] += dZj;
                    Y[j][celli] -= dZj;
                }
            }
        }
    }

    // Condensation

    if (condensation_->modelType() != "none")
    {
        PtrList<scalarField> rhoCont(thermo_.rhoCont(contSpecies));

        forAll(rho, celli)
        {
            const conData cdata
            (
                condensation_->rate
                (
                    dcm[celli],
                    p[celli],
                    T[celli],
                    entryList(Y,celli),
                    entryList(Z,celli),
                    entryList(pSat,celli),
                    entryList(gamma,celli),
                    entryList(D,celli),
                    entryList(rhoCont,celli),
                    entryList(rhoDisp,celli),
                    entryList(sigma,celli)
                )
            );

            if (cdata.active())
            {
                scalarList M0(dist.size(), 0.0);

                scalar sumM(0.0);

                forAll(sections, i)
                {
                    M0[i] = max(sections[i].M().field()[celli],0.0);

                    sumM += M0[i];

                    sections[i].M().field()[celli] = 0.0;
                }

                const scalarList Y0(entryList(Y,celli));
                const scalarList Z0(entryList(Z,celli));

                const scalar d(min(dcm[celli],dMax_));

                scalar dAlpha(0.0);

                forAll(activeSpecies, j)
                {
                    const scalar a((Y0[j]+Z0[j])*cdata.source()[j]);

                    const scalar b
                    (
                        max
                        (
                            (cdata.source()[j]+cdata.sink()[j]),
                            VSMALL
                        )
                    );

                    Z[j][celli] = max
                    (
                        a/b
                      + (Z0[j] - a/b)
                        * Foam::exp(-b*sumM*d/rDeltaT[celli]),
                        0.0
                    );

                    Y[j][celli] = max(Y0[j]+Z0[j]-Z[j][celli],0.0);

                    dAlpha += (Z[j][celli]-Z0[j]);

                    I_[j].field()[celli] =
                        rho[celli]*(Z[j][celli]-Z0[j])*rDeltaT[celli];
                }

                // I(s) ~ const

                const scalar Gamma
                (
                    dAlpha*rDeltaT[celli]
                  / (max(d*sumM,VSMALL))
                );

/*
                // I(s)/d ~ const

                const scalar Gamma
                (
                    dAlpha
                  * Foam::pow(6.0/(rhol[celli]*pi), 1.0/3.0)*rDeltaT[celli]
                  / (max(d*sumM,VSMALL))
                );
*/

                forAll(sections, i)
                {
                    if (M0[i] > SMALL)
                    {
                        // I(s) ~ const

                        // Note: s can become negative

                        const scalar s
                        (
                            dist[i].x()
                            + Gamma/rDeltaT[celli]*dist[i].d(rhol[celli])
                        );

/*
                        // I(s)/d ~ const

                        // Note: without clipping the argument of the 3/2 power
                        // can become negative. The clipping causes an
                        // inconsistency in the mass transfer in f

                        const scalar s
                        (
                            Foam::pow
                            (
                                max
                                (
                                    Foam::pow(dist[i].x(), 2.0/3.0)
                                    + 2.0/3.0*Gamma/rDeltaT[celli],
                                    0.0
                                ),
                                3.0/2.0
                            )
                        );
*/

                        if (s >= dist.xMin())
                        {
                            interp.addToM(s, M0[i], celli);
                        }
                    }
                }
            }
        }
    }

// Coalescence

if (coalescence_->modelType() != "none")
{
    const scalar gMag = max(mag(g_), SMALL);

    // Detect turbulence on/off
    const IOdictionary& turbProps =
        db().lookupObject<IOdictionary>("turbulenceProperties");

    const word simulationType(turbProps.get<word>("simulationType"));
    const bool turbulenceOn = (simulationType == "RAS" || simulationType == "LES");

    // --------------------------------------------------
    // epsilon estimation
    // --------------------------------------------------

    const volScalarField* epsilon = nullptr;

    static autoPtr<volScalarField> epsilonTmpPtr;      // for kOmega/SST
    static autoPtr<volScalarField> epsilonDefaultPtr; // fallback

    if (turbulenceOn)
    {
        // 1) k-epsilon, LES, or user-provided epsilon
        if (db().foundObject<volScalarField>("epsilon"))
        {
            epsilon = &db().lookupObject<volScalarField>("epsilon");
        }
        // 2) k-omega / k-omega-SST: build epsilon = betaStar * k * omega
        else if
        (
            db().foundObject<volScalarField>("omega") &&
            db().foundObject<volScalarField>("k")
        )
        {
            const volScalarField& omega =
                db().lookupObject<volScalarField>("omega");

            const volScalarField& k =
                db().lookupObject<volScalarField>("k");

            const scalar betaStar = 0.09;

            if (!epsilonTmpPtr.valid())
            {
                epsilonTmpPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "epsilonFromOmega",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        betaStar * omega * k
                    )
                );
            }
            else
            {
                epsilonTmpPtr() = betaStar * omega * k;
            }

            epsilon = epsilonTmpPtr.get();
        }
    }

    // 3) fallback epsilon (missing turbulence fields)
    if (!epsilon)
    {
        if (!epsilonDefaultPtr.valid())
        {
            epsilonDefaultPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilonDefault",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "epsilon",
                        dimensionSet(0, 2, -3, 0, 0, 0, 0),
                        1e-6
                    )
                )
            );
        }

        epsilon = epsilonDefaultPtr.get();

        if (turbulenceOn)
        {
        WarningInFunction
            << "epsilon not found in turbulence model. "
            << "Using default epsilon = 1e-6" << endl;
        }
    }

    // --------------------------------------------------
  
        const scalarField mug(thermo_.thermoCont().mu());
        const scalarField rhog(thermo_.thermoCont().rho());
        
        const scalar pi = constant::mathematical::pi;
        const scalar kB = constant::physicoChemical::k.value();
        const scalar NA = constant::physicoChemical::NA.value();
               
      
    	
    	const scalar R = constant::physicoChemical::R.value();

        const PtrList<coalescencePair>& pairs =
            system_->coalescencePairs();

        if (pairs.size() == 0)
        {
            system_->generateCoalescencePairs();
        }
        
         // Choose particle shape type: "fractal"/"custom" or "solid"
	       const word& particleShape = particleShape_;

        forAll(rho, celli)
        {
                scalarList d(sections.size(), 0.0);
                scalarList M0(dist.size(), 0.0);
                
                scalarList N(sections.size(), 0.0); // Number of monomers per aggregate
                scalarList volEqDiam(sections.size(), 0.0);  // Volume-equivalent diameter
                
                scalarList ds_list(sections.size(), 0.0);
                scalarList Df_list(sections.size(), 0.0);
                
                forAll(sections, i)
                {
                    d[i] = dist[i].d(rhol[celli]);
                    M0[i] = max(sections[i].M().field()[celli], 0.0);
                    
				if (particleShape == "fractal" || particleShape == "custom")
					{
					// -------------------------------
					// Fractal/aggregate calculations
					// -------------------------------
							
					const scalar ds = 2.0 * monoRad();
					const scalar Df = frDim();
					const scalar mono_rho = monoRho();
							
					scalar V_agg     = (pi / 6.0) * pow(d[i], 3);
					scalar V_monomer = (pi / 6.0) * pow(ds, 3);

					N[i] = (V_agg / V_monomer) * (rhol[celli] / mono_rho);   //number of monomers N_i = M_agg / M_monomer
					volEqDiam[i] = pow(6.0 * V_agg / pi, 1.0 / 3.0);
							
					ds_list[i] = ds;
					Df_list[i] = Df;
							
				
					}
				else if (particleShape == "solid")
					{
					// -------------------------------
					// Pure spherical particle
					// -------------------------------
					N[i] = 1.0;  // not used, but set to avoid uninitialized values
					volEqDiam[i] = d[i];  // volume-equivalent diameter = actual diameter
							
					ds_list[i] = 0.0;   // dummy
                                       Df_list[i] = 0.0;   // dummy
					}
		} 

            
                forAll(pairs, k)
                {
                    const coalescencePair& pair = pairs[k];

                    const label i(pair.i());
                    const label j(pair.j());
                    
                  
                  // ---------------------------
    	          // Diameters & Mass
                 // ---------------------------
			scalar dfi = VSMALL, dfj = VSMALL;   // fractal diameters
			scalar dmi = VSMALL, dmj = VSMALL;   // mobility diameters
			scalar mi  = VSMALL, mj  = VSMALL;   // particle masses
                      
                    // ==============================================
                		// Fractal/Custom And Solid Diameter Calculations (Jacobson Model)
                   // ==============================================


    		if (particleShape == "fractal" || particleShape == "custom")
    			{
        			// === FRACTAL MODEL ===
        		dfi = computeFractalDiameter(ds_list[i], N[i], Df_list[i]);
        		dfj = computeFractalDiameter(ds_list[j], N[j], Df_list[j]);

        		dmi = computeMobilityDiameter(ds_list[i], N[i], Df_list[i], dfi);
        		dmj = computeMobilityDiameter(ds_list[j], N[j], Df_list[j], dfj);

        		mi = (pi / 6.0) * rhol[celli] * pow(ds_list[i], 3) * N[i]; // mass of fractal = N * mass of monomar
        		mj = (pi / 6.0) * rhol[celli] * pow(ds_list[j], 3) * N[j];
    			}
    		else if (particleShape == "solid")
    			{
        			// === PURE SPHERICAL MODEL ===
        		dfi = d[i];
        		dfj = d[j];

        		dmi = d[i];
        		dmj = d[j];

        		mi = (pi / 6.0) * rhol[celli] * pow(d[i], 3);
        		mj = (pi / 6.0) * rhol[celli] * pow(d[j], 3);
    			}

                    
                    // ==============================================
                		// Fractal Diameter Calculations (Jacobson Model)
                   // ==============================================


             scalar m_gas = rhog[celli] * R * T[celli] / p[celli] / NA;

            // Air mean free path
            
            scalar lambda = mug[celli] / p[celli] * sqrt(pi * kB * T[celli] / (2.0 * m_gas));
            
            // Knudsen number for fractals

            scalar kni = 2.0 * lambda / dmi;
            scalar knj = 2.0 * lambda / dmj;
            
            // Cunningham slip correction for fractals
            
            scalar Cci = 1.0 + (kni/2.0) * (2.34 + 1.05 * exp(-0.39 / (kni/2.0)));
            scalar Ccj = 1.0 + (knj/2.0) * (2.34 + 1.05 * exp(-0.39 / (knj/2.0)));
            
            // Particle diffusion coefficient for fractals
            
            scalar Di = kB * T[celli] * Cci / (3.0 * pi * mug[celli] * dmi);
            scalar Dj = kB * T[celli] * Ccj / (3.0 * pi * mug[celli] * dmj);
            
            // Thermal speeds of fractals
            
            scalar vpi = sqrt(8.0 * kB * T[celli] / (pi * mi));
            scalar vpj = sqrt(8.0 * kB * T[celli] / (pi * mj));
            
            scalar va = mug[celli]/rhog[celli];  // kinematic viscosity of fluid
            
            // Transition parameters (using dmi/dmj)
            
            // fractal mean free path

            scalar lambdapi = 8.0 * Di / (pi * vpi);
            scalar lambdapj = 8.0 * Dj / (pi * vpj);
            
            // mean distance calculation 
            
            scalar deltai = (pow(dmi + lambdapi, 3) - pow(sqr(dmi) + sqr(lambdapi), 1.5)) / (3.0 * dmi * lambdapi) - dmi;
            scalar deltaj = (pow(dmj + lambdapj, 3) - pow(sqr(dmj) + sqr(lambdapj), 1.5)) / (3.0 * dmj * lambdapj) - dmj;
            
            // Denominator terms (using dfi/dfj for collision cross-section)

            scalar dterm1 = (dfi + dfj) / (dfi + dfj + 2.0 * sqrt(sqr(deltai) + sqr(deltaj))); // first term in denominator
            scalar dterm2 = 8.0 * (Di + Dj) / ((dfi + dfj) * sqrt(sqr(vpi) + sqr(vpj))); // second term in denominator
            
            // Final Fuchs kernel for fractals
            
            scalar beta_Fuchs = 2.0 * pi * (Di + Dj) * (dfi + dfj) / (dterm1 + dterm2); // KB term
            
            scalar beta = beta_Fuchs; // Default: pure Brownian Fuchs kernel (for laminar)
            
            scalar beta_turb_inertial = 0.0;
            scalar beta_turb_shear = 0.0;
            
            if (turbulenceOn)
            {
                const scalar epsilonVal = epsilon->internalField()[celli];
                
                
                
                if (epsilonVal > SMALL)
                {
                    
                    scalar vfi = (rhol[celli] - rhog[celli]) * sqr(dmi) *gMag* Cci / (18.0 * mug[celli]); // Terminal fall speed
                    scalar vfj = (rhol[celli] - rhog[celli]) * sqr(dmj) *gMag* Ccj / (18.0 * mug[celli]); // Terminal fall speed

                    beta_turb_inertial = pi * pow(epsilonVal, 0.75) * sqr(dfi + dfj) * fabs(vfi - vfj) / (4.0 *gMag* pow(va, 0.25));
                    
                    beta_turb_shear = 0.125 * sqrt(8.0 * pi * epsilonVal / (15.0 * va)) * pow(dfi + dfj, 3);
                    
                    //total beta
                    
                    beta += beta_turb_inertial + beta_turb_shear;
                }
            }
/*                   
                   betaOutput 
    			<< celli << " "
   			<< i << " "
    			<< j << " "
    			<< volEqDiam[i] << " "
    			<< volEqDiam[j] << " "
    			<< beta_Fuchs << " "
    			<< beta_turb_inertial << " "
    			<< beta_turb_shear << " "
    			<< beta << endl;

*/
                    scalar& Mi = sections[i].M().field()[celli];
                    scalar& Mj = sections[j].M().field()[celli];

                    const scalar f
                    (
                        min
                        (
                            M0[i]*M0[j]*rho[celli]*beta/rDeltaT[celli],
                            min(Mi, Mj)
                        )
                    );

                    Mi -= f;
                    Mj -= f;

                    interp.addToM(pair.idata(), pair.s(), f, celli);
                }
            }
        
    }

    system_->rescale();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModels::fixedSectional::fixedSectional
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    aerosolModel(modelType, mesh, aerosolProperties),
    
    g_
    (
        uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedVector("g", dimVelocity/dimTime, vector::zero)
        ).value()
    ),
    system_(),
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("J", dimless/dimVolume/dimTime, 0)
    ),
    I_(thermo_.activeSpecies().size()),
    rescale_(coeffs_.lookupOrDefault<Switch>("rescale", true))
    

{
    system_.reset(new fixedSectionalSystem(*this, coeffs()));

    const speciesTable& activeSpecies = thermo_.activeSpecies();

    forAll(activeSpecies, j)
    {
        I_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word(activeSpecies[j]+":I"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("I", dimMass/dimTime/dimVolume, 0)
            )
        );
    }

    phiEff_.setSize(system_->distribution().size());

    forAll(system_->distribution(), i)
    {
        const section& sec = system_->distribution()[i];

        drift_->fields().add(sec.M());

        mesh.setFluxRequired(sec.M().name());

        phiEff_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("phiEff", sec.M().name()),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("phi", dimless/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModels::fixedSectional::~fixedSectional()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aerosolModels::fixedSectional::correct()
{
    // Only the M-field must be initialized at startup. Sources are corrected in
    // the solvePost step directly

    forAll(system_->distribution(), i)
    {
        section& sec = system_->distribution()[i];

        if(!sec.validM())
        {
            sec.initM(word(coeffs().lookup("initFromPatch")));
        }

        // needed for restart 
        sec.M().correctBoundaryConditions();
    }
}

void Foam::aerosolModels::fixedSectional::correctDriftFlux()
{
    drift_->correct();
    updateDrift();
}


void Foam::aerosolModels::fixedSectional::solvePre()
{}


void Foam::aerosolModels::fixedSectional::solvePost()
{
    solveSpatial();
    solveInternal();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::aerosolModels::fixedSectional::R(const volScalarField& Y) const
{
    tmp<volScalarField> I
    (
        new volScalarField
        (
            IOobject
            (
                Y.name() + ":I:zero",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("I", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    return fvm::Su(I, Y);
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::Qdot() const
{
    return condensation_->Qdot(I_);
}


Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::meanDiameter
(
    const scalar p,
    const scalar q
) const
{
    dimensionedScalar dMin("d", dimLength, dMin_);
    dimensionedScalar dMax("d", dimLength, dMax_);

    return max(min(system_->meanDiameter(p,q),dMax),dMin);
}

Foam::tmp<Foam::scalarField>
Foam::aerosolModels::fixedSectional::meanDiameter
(
    const scalar p,
    const scalar q,
    const label patchi
) const
{
    return max(min(system_->meanDiameter(p,q,patchi),dMax_),dMin_);
}

Foam::tmp<Foam::volScalarField>
Foam::aerosolModels::fixedSectional::medianDiameter
(
    const scalar p
) const
{
    dimensionedScalar dMin("d", dimLength, dMin_);
    dimensionedScalar dMax("d", dimLength, dMax_);

    return max(min(system_->medianDiameter(p),dMax),dMin);
}

void Foam::aerosolModels::fixedSectional::clearRates()
{
    J_ *= 0.0;

    forAll(I_, j)
    {
        I_[j] *= 0.0;
    }
}

bool Foam::aerosolModels::fixedSectional::read()
{
    if (aerosolModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
