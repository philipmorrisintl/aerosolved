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

{
    // Prepare

    volScalarField Yt(0.0*Y[0]);

    const volScalarField mut(turbulence->mut());

    const surfaceScalarField& phiDrift = aerosol->phiDrift();
    const surfaceScalarField& phiInertial = aerosol->phiInertial();
    const surfaceScalarField& phiBrownian = aerosol->phiBrownian();

    const fv::convectionScheme<scalar>& mvPhi = aerosol->mvPhi();
    const fv::convectionScheme<scalar>& mvPhiDrift = aerosol->mvPhiDrift();
    const fv::convectionScheme<scalar>& mvPhiInertial = aerosol->mvPhiInertial();
    const fv::convectionScheme<scalar>& mvPhiBrownian = aerosol->mvPhiBrownian();

    phiEffY[inertIndex] = phi;

    // Continuous mass fraction equations

    forAll(Y, j)
    {
        if (j != inertIndex)
        {
            volScalarField& Yj = Y[j];

            const volScalarField& D = aerosol->DCont()[j];

            fvScalarMatrix YjEqn
            (
                fvm::ddt(rho, Yj)
              + mvPhi.fvmDiv(phi, Yj)
              - mvPhiDrift.fvmDiv(phiDrift, Yj)
             ==
                aerosol->R(Yj)
              + fvm::laplacian(mut+rho*D, Yj,"laplacian(mut+rho*D,Yj)")
              + fvOptions(rho, Yj)
            );

            YjEqn.relax();

            fvOptions.constrain(YjEqn);

            YjEqn.solve(mesh.solver("Y"));

            fvOptions.correct(Yj);

            Yj.max(0.0);
            Yt += Yj;

            phiEffY[j] = YjEqn.flux();

            phiEffY[inertIndex] -= phiEffY[j];
        }
    }

    // Dispersed mass fraction equations

    const volScalarField& D = aerosol->DDisp();

    forAll(Z, j)
    {
        volScalarField& Zj = Z[j];

        fvScalarMatrix ZjEqn
        (
            fvm::ddt(rho, Zj)
          + mvPhi.fvmDiv(phi, Zj)
          + mvPhiInertial.fvmDiv(phiInertial, Zj)
          + mvPhiBrownian.fvmDiv(phiBrownian, Zj)
          - mvPhiDrift.fvmDiv(phiDrift, Zj)
          ==
            aerosol->R(Zj)
          + fvm::laplacian(mut+rho*D, Zj,"laplacian(mut+rho*D,Zj)")
          + fvOptions(rho, Zj)
        );

        ZjEqn.relax();

        fvOptions.constrain(ZjEqn);

        ZjEqn.solve(mesh.solver("Y"));

        fvOptions.correct(Zj);

        Zj.max(0.0);
        Yt += Zj;

        phiEffZ[j] = ZjEqn.flux();

        phiEffY[inertIndex] -= phiEffZ[j];
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    Y[inertIndex].correctBoundaryConditions();
}
