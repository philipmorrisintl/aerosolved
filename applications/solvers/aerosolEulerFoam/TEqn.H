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
    const surfaceScalarField phiv(fvc::absolute(phi/fvc::interpolate(rho), U));

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(phi, T)
      ==
        fvm::laplacian(thermo.alphaEff(turbulence->mut()), T,"laplacian(alphaEff,T)")
      - (
            fvc::div(phiv, p, "div(phiv,p)")
          + fvc::ddt(rho, K)
          + fvc::div(phi, K)
          - aerosol->Qdot()
        )
      / thermo.Cv()
      + fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    thermo.correctThermo();
    thermo.correct();

    Info<< "min(T) = " << min(T).value() << ", max(T) = " << max(T).value() << endl;
}
