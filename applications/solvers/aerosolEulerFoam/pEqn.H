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

\file pEqn.H
\brief Solve the \f$p\f$-equation.

In this file a number of things are done:

- The pressure equation is formulated. For the \f$k\f$th iteration (or \f$(k-1)\f$th
  corrector) the pressure equation is given by Eq. (3.34) in the thesis. Depending
  on the fluidThermo model which is selected, the ddt-term may be present implicitly,
  explicitly, or not present at all in case of incompressible flow.

- The pressue equation is solved. For the final non-orthogonal iteration the flux
  field \f$\phi\f$ is updated from the pressure equation. This update corresponds to
  Eq. (3.41).

- The velocity \f$\mathbf{u}\f$ is updated. This corresponds to Eq. (3.40).

- The material time derivative of pressure \f$\mathrm{D}_tp\f$ is updated for further
  use in the temperature equation.

*/

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(rho*HbyA) & mesh.Sf())
);

while (piso.correctNonOrthogonal())
{
    fvScalarMatrix p1Eqn
    (
        fvc::div(phiHbyA + phidRho)
      - fvm::laplacian(rhorAUf, p1)
    );

    if (thermo.type() == "semiIncompressible")
    {
        p1Eqn = (p1Eqn + fvc::ddt(rho));
    }
    else if (thermo.type() == "compressible")
    {
        p1Eqn = (p1Eqn + fvm::ddt(psi, p1) + p0 * fvc::ddt(psi));
    }
    else if (thermo.type() == "incompressible")
    {
        // No ddt(rho)
    }

    p1Eqn.solve();

    // Update the initial residual if we're in the first non-orthogonal corrector

    if (piso.corrNonOrtho() == 1)
    {
        piso.updateResidual();
    }

    // Update the flux in the final non-orthogonal corrector, Eq. (3.32) or (3.41)

    if (piso.finalNonOrthogonalIter())
    {
        phi == phiHbyA + p1Eqn.flux();
    }
}

// Update velocity, Eq. (3.31) or (3.40)

U = rAU*(UEqn.H() - fvc::grad(p1));

U.correctBoundaryConditions();

DpDt == fvc::ddt(p1)+(U & fvc::grad(p1));
