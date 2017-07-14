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

\file MEqn.H
\brief Solve the \f$M_i\f$-equation's first fractional step

In this file the \f$M_i\f$-equation's first fractional step is solved. This
corresponds to Eq. (3.25) for the implicit predictor step and Eq. (3.44)
for the corrector steps (see \cite thesis) for the set \f$\mathbf{X}={M_i}\f$. The
Crank-Nicolson scheme is implemented manually. This means that the ddt-term
reduces to the default Euler scheme (adding a \f$1/\Delta t\f$ to the
diagonal and \f$M_i^{m}/\Delta t\f$ to the source. The fluxes are always
taken explicitly and are based on fluxes.H.

*/

if (piso.corr() == 0 && piso.theta() < 1.0)
{
    forAll(M, i)
    {
        SddtM[i] = fvc::div(phiM[i]);
    }
}

forAll(M, i)
{
    volScalarField& Mi = M[i];
    const volScalarField& SddtMi = SddtM[i];
    const volScalarField& Ji = J[i];

    fvScalarMatrix MiEqn
    (
        fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rhoStar,Mi)
      + im * fvc::div(phiM[i])
      + ex * SddtMi
      ==
        Ji
    );

    Mi = MiEqn.H() / MiEqn.A();
    Mi.correctBoundaryConditions();
}
