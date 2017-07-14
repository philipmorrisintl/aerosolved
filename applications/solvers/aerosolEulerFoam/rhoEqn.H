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

\file rhoEqn.H
\brief Solve the \f$\rho\f$-equation

For the first iteration an explicit extrapolation is performed of the continuity
equation Eq. (3.22) in \cite thesis. Also a new \f$\rho\f$-equation is formulated
but now for \f$\rho^\star\f$. This alternative density 'feels' the updated flux
of the phidRho. \f$\rho^\star\f$ should be very close to \f$\rho\f$, but ensures
that mass fractions will always sum up to unity up to machine precision. This is
beneficial for numerical consistency.

*/

if (piso.corr() == 0)
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi + phidRho)
    );

    rhoEqn.solve(mesh.solver("rho"));

    // Re-evaluated fluxes for explicit CN source term

    if (piso.theta() < 1.0)
    {
        SddtRhoStar = - fvc::div(phi + phidRho);
    }

    rhoStar == rho;
    rhoStar.oldTime() == rho.oldTime();
}

fvScalarMatrix rhoStarEqn
(
    fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rhoStar)
  + im * fvc::div(phi + phidRho)
  ==
  + ex * SddtRhoStar
);

rhoStar = rhoStarEqn.H() / rhoStarEqn.A();
rhoStar.correctBoundaryConditions();
