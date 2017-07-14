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

/**

\file aerosolEulerFoam.C
\brief Main aerosol PISO solver

This is the main source file of the aerosolEulerFoam solver. It is based
on \cite thesis, in particular Chptr. 3 in that work, which, in turn, is
mainly based on Issa's compressible PISO algorithm for reacting flows
\cite Issa:1991aa.

The executable has four options:

- -timeAveraing: enable time averaging.

- -plausibilityCheck: checks if the solution is still plausible as
  defined by the plausibilityLimits dictionary.

- -externalGradP: reads the external pressure gradient from the
  transportProperties dictionary and inserts it into the
  momentum equation.

- -positivity: enables clipping to zero to ensure positivity of the
  solution.

References to equation numbers in comments in the code are with respect to
\cite thesis.

*/

// OpenFOAM includes

#include "multivariateScheme.H"
#include "fvCFD.H"
#include "EulerDdtScheme.H"

// AeroSolved includes

#include "aerosolModel.H"
#include "fluidThermo.H"
#include "driftVelocityModel.H"
#include "nucleationModel.H"
#include "condensationEvaporationModel.H"
#include "diffusionModel.H"
#include "brownianDiffusionModel.H"
#include "conductivityModel.H"
#include "coalescenceModel.H"
#include "viscosityModel.H"
#include "turbModel.H"
#include "pisoControl.H"
#include "timeAveraging.H"
#include "plausibility.H"

int main(int argc, char *argv[])
{
    argList::validOptions.insert("timeAveraging","");
    argList::validOptions.insert("plausibilityCheck","");
    argList::validOptions.insert("externalGradP","");
    argList::validOptions.insert("positivity","");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    pisoControl piso(mesh);

    // Create the plausibility object

    Plausibility plausibility
    (
        thermo,
        runTime,
        args.options().found("plausibilityCheck")
    );

    // Get theta for the custom theta scheme implementations in Yj, Zj and Mi.

    const scalar theta = piso.theta();
    scalar im = theta;
    scalar ex = (1.0-theta);

    Info << "Using theta = " << theta << " (implicit = "
         << im << ", explicit = " << ex << ")" << endl;

    // Print the selected schemes

    Info << endl << mesh.schemesDict() << endl << endl;

    Info<< "\nStarting time loop\n" << endl;

    // Start the time loop

    while (runTime.run())
    {
        #include "readTimeControls.H"

        // Compute CFL numbers and restrict the time step

        #include "stability.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ---- Predictor step  ----------------------------------------------

        if (doEuler)
        {
            im = 1.0;
            ex = 0.0;

            doEuler = false;
        }
        else
        {
            im = theta;
            ex = (1.0-theta);
        }

        // Updated explicit fluxes phiM, phiY, phiZ and phidRho

        #include "fluxes.H"

        // Update the aerosol model

        aerosol.update();

        // Explicit prediction of density, Eq. (3.22)

        #include "rhoEqn.H"

        // Implicit prediction of temperature, Eq. (3.24)

        #include "TEqn.H"

        // First fractional step of the set X, Eq. (3.25)

        #include "YEqn.H"
        #include "ZEqn.H"
        #include "MEqn.H"

        // Second fractional step of the set X, Eq (3.26)

        aerosol.fractionalStepInternal();

        // Update psi based on the first predictions

        thermo.updatePsi();

        #include "coeffs.H"

        // First implicit prediction of the velocity, Eq. (3.23)

        #include "UEqn.H"

        // ---- Corrector steps  ---------------------------------------------

        while(piso.loop())
        {
            if (piso.corr() > 1)
            {
                // Update psi and density

                thermo.updatePsi();
                thermo.updateRho();

                #include "coeffs.H"
            }

            // Pressure equation, Eq. (3.420

            #include "pEqn.H"

            // Update fluxes phiM, phiY, phiZ and phidRho

            #include "fluxes.H"

            thermo.updateRho();

            #include "rhoEqn.H"

            // Correct the temperature, Eq. (3.43)

            #include "TEqn.H"

            // First fractional step of the set X, Eq. (3.44)

            #include "YEqn.H"
            #include "ZEqn.H"
            #include "MEqn.H"

            // Second fractional step of the set X, Eq. (3.45)

            aerosol.fractionalStepInternal();
        }

        // Third fractional step (coagulation)

        aerosol.fractionalStepExternal();

        // Rescale solution

        if(args.options().found("positivity"))
        {
            thermo.rescale(Y, Z, M, true, true);
        }
        else
        {
            thermo.rescale(Y, Z, M, true, false);
        }

        // Correct the size distribution to satisfy the consistency relation, Eq. (2.23) or (2.31)

        aerosol.correctSizeDistribution();

        aerosol.checkConsistency();

        // Time averaging

        if(args.options().found("timeAveraging"))
        {
            timeAveraging
            (
                meanDataFields, U, thermo, aerosol, startAveraging, runTime
            );
        }

        // Plausibility check

        if (args.options().found("plausibilityCheck"))
        {
            plausibility.check();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    return 0;
}

// ************************************************************************* //
