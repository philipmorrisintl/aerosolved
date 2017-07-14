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

\file ZEqn.H
\brief Solve the \f$Z_j\f$-equations

For the first iteration, update the explicit source term

For all species \f$j\f$, formulate the \f$Z_j\f$-equation and solve. For the
\f$k\f$th iteration this corresponds to Eq. (3.44) with the set
\f$\mathbf{X}={Z_j}\f$. The time integration is hard-coded to be the
\f$\theta\f$-scheme, so that the ddt-term must me set using the Euler scheme.

*/

if (piso.corr() == 0 && piso.theta() < 1.0)
{
    forAll(species, j)
    {
        SddtZ[j] = fvc::div(phiZ[j]);
    }
}

forAll(species, j)
{
    volScalarField& Zj = Z[j];
    const volScalarField& Sj = S[j];
    volScalarField& SddtZj = SddtZ[j];

    fvScalarMatrix ZjEqn
    (
        fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rhoStar,Zj)
      + im * fvc::div(phiZ[j])
      + ex * SddtZj
      ==
        Sj
    );

    Zj = ZjEqn.H() / ZjEqn.A();
    Zj.correctBoundaryConditions();
}
