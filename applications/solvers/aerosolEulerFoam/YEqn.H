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

\file YEqn.H
\brief Solve the \f$Y_j\f$-equations

For all species \f$j\f$, formulate the \f$Y_j\f$-equation and solve. For the
\f$k\f$th iteration this corresponds to Eq. (3.44) with the set
\f$\mathbf{X}={Y_j}\f$. The time integration is hard-coded to be the
\f$\theta\f$-scheme, so that the ddt-term must me set using the Euler scheme.

*/

if (piso.corr() == 0)
{
    diffusion.update();

    forAll(species, j)
    {
        SddtYD[j] == fvc::laplacian(rho*DY[j], Y[j])
                   - fvc::div(phic, Y[j], "div(phic,Y)");

        if (piso.theta() < 1.0)
        {
            SddtY[j] = fvc::div(phiY[j]);
        }
    }
}

forAll(species, j)
{
    volScalarField& Yj = Y[j];
    const volScalarField& Sj = S[j];
    volScalarField& SddtYj = SddtY[j];

    fvScalarMatrix YjEqn
    (
        fv::EulerDdtScheme<scalar>(mesh).fvmDdt(rhoStar,Yj)
      + im * fvc::div(phiY[j])
      + ex * SddtYj
      ==
      - Sj
    );

    YjEqn = (YjEqn == SddtYD[j]);

    Yj = YjEqn.H() / YjEqn.A();
    Yj.correctBoundaryConditions();
}
