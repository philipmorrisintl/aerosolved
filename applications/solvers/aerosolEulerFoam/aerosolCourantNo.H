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

const scalar maxAerosolCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxAerosolCo", maxCo);

const scalar maxDriftCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxDriftCo", maxCo);

scalar aerosolCoNum = 0.0;
scalar driftCoNum = 0.0;
scalar relCoNum = 0.0;

if (aerosol->drift().inertial().type() != "none")
{
    forAll(Z, j)
    {
        aerosolCoNum = max
        (
            gMax
            (
                mag
                (
                    aerosol->R(Z[j])().source()
                  / (rho.field()*mesh.V().field())
                  * runTime.deltaTValue()
                )
            ),
            aerosolCoNum
        );
    }

    const surfaceScalarField& phiInertial = aerosol->phiInertial();
    const surfaceScalarField& phiBrownian = aerosol->phiBrownian();

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi+phiInertial+phiBrownian))().primitiveField()
      / rho.primitiveField()
    );

    driftCoNum =  0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    scalarField sumPhiv
    (
        fvc::surfaceSum(mag(phiInertial+phiBrownian))().primitiveField()
      / rho.primitiveField()
    );

    relCoNum = 0.5*gMax(sumPhiv/mesh.V().field())*runTime.deltaTValue();
}

Info<< "Aerosol Courant Number max: " << aerosolCoNum << endl;
Info<< "Drift Courant Number max: " << driftCoNum << endl;
Info<< "Relative Courant Number max: " << relCoNum << endl;

// ************************************************************************* //
