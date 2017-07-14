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

\file TEqn.H
\brief Solve the \f$T\f$-equation

For the specified number of non-orthogonal correctors the temperature equation
is formulated and solved. The time discretization is taken from the fvSchemes
file. For the \f$k\f$th iteration the formulated equation corresponds to Eq.
(3.43) from \cite thesis.

*/

while(piso.correctNonOrthogonal())
{
    fvScalarMatrix TEqn
    (
        CpEff*fvm::ddt(rho, T)
      + CpEff*fvm::div(phi, T, "div(phi,T)")
      ==
        fvm::laplacian(kEff+kTurb, T)
      + (muEff+muTurb) * (fvc::grad(U) && fvc::grad(U))
      + DpDt
    );

    TEqn.solve();
}
