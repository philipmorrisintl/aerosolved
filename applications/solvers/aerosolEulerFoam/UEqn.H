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

\file UEqn.H
\brief Solve the \f$\mathbf{u}\f$-equation

Formulate the \f$\mathbf{u}\f$-equation (see Eq. (3.39) in \cite thesis for the
\f$k\f$th corrector) and solve it using the old pressure in the first iteration
(i.e., the implicit predictor step of the PISO algorithm). For subsequent
iterations the UEqn object must be formulated as it is required for the pEqn
file. An external pressure gradient may also be added.

*/

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  ==
    fvm::laplacian(muEff+muTurb, U)
  + fvc::div((muEff+muTurb)*dev2(Foam::T(fvc::grad(U))), "div(mu*dev2(T(grad(U))))")
  + fvc::grad(turbulence.K(U))
  + externalGradP
);

if (piso.corr() == 0)
{
    solve(UEqn == - fvc::grad(p1));
}
