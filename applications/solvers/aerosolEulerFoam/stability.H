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

\file stability.H
\brief Compute and show all stability numbers

The following stability numbers are computed (mean and max):

- Drift CFL number, taking the maximum CFL number of all sectional
  \f$\mathbf{v}_i\f$-fields.
- Droplet diffusional CFL number, taking the maximumum CFL-like number of all
  the sectional diffusional velocities.
- Vapor diffusional CFL number, taking the maximum CFL-like number of all the
  vapor diffusional velocities.
- Vapor diffusion correction flux (i.e., \f$\mathbf{u}_c\f$ in Eq. (2.19) in
  \cite thesis), taking the CFL number of this correction flux.
- Viscous CFL number based on the viscous flux.
- Standard CFL number based on the flux field \f$\phi\f$.

All these numbers are printed to the screen.

*/

scalar driftCoNum = 0.0;
scalar driftMeanCoNum = 0.0;

scalar diffMCoNum = 0.0;
scalar diffMMeanCoNum = 0.0;

scalar diffYCoNum = 0.0;
scalar diffYMeanCoNum = 0.0;

scalar diffCorrYCoNum = 0.0;
scalar diffCorrYMeanCoNum = 0.0;

scalar muCoNum = 0.0;
scalar muMeanCoNum = 0.0;

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

CFL = dimensionedScalar("zero", dimless, 0.0);

if (mesh.nInternalFaces())
{
    // Mixture velocity

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().internalField()
      / rho.internalField()
    );

    CFL.internalField() = 0.5*sumPhi/mesh.V().field()*runTime.deltaTValue();

    CoNum = gMax(CFL.internalField());

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    // Vapor diffusion correction flux

    scalarField sumPhic
    (
        fvc::surfaceSum(mag(phic))().internalField()
      / rho.internalField()
    );

    diffCorrYCoNum = gMax(0.5*sumPhic/mesh.V().field()*runTime.deltaTValue());

    diffCorrYMeanCoNum =
        0.5*(gSum(sumPhic)/gSum(mesh.V().field()))*runTime.deltaTValue();

    // Droplet velocity

    forAll (phid, i)
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi+phid[i]))().internalField()
          / rho.internalField()
        );

        scalarField CFLi(0.5*sumPhi/mesh.V().field()*runTime.deltaTValue());

        forAll(CFL.internalField(), j)
        {
            CFL.internalField()[j] = max(CFL.internalField()[j], CFLi[j]);
        }

        driftCoNum = max(gMax(CFLi), driftCoNum);

        driftMeanCoNum +=
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    driftMeanCoNum /= scalar(max(aerosol.P(), 1));

    scalarField sumSurf
    (
        fvc::surfaceSum(mesh.magSf()*mesh.surfaceInterpolation::deltaCoeffs())().internalField()
    );

    // Droplet diffusion

    diffMCoNum = gMax(DM[0].internalField()*sumSurf/mesh.V().field())*runTime.deltaTValue();

    diffMMeanCoNum = gSum(DM[0].internalField()*sumSurf)/gSum(mesh.V().field())*runTime.deltaTValue();

    // Vapor diffusion

    forAll(species, j)
    {
        diffYCoNum = max
        (
            diffYCoNum,
            gMax(DY[j].internalField()*sumSurf/mesh.V().field())*runTime.deltaTValue()
        );

        diffYMeanCoNum = max
        (
            diffYMeanCoNum,
            gSum(DY[j].internalField()*sumSurf)/gSum(mesh.V().field())*runTime.deltaTValue()
        );
    }

    // Viscosity

    muCoNum = gMax(muEff.internalField()*sumSurf/mesh.V().field())*runTime.deltaTValue();

    muMeanCoNum = (gSum(muEff.internalField()*sumSurf)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info << "Restricting stability numbers:" << endl << endl;

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

Info<< "Drift Courant Number mean: " << driftMeanCoNum
    << " max: " << driftCoNum << endl;

Info<< "Droplet Diffusion Stability Number mean: " << diffMMeanCoNum
    << " max: " << diffMCoNum << endl;

Info<< "Vapor Diffusion Correction Stability Number mean: " << diffCorrYMeanCoNum
    << " max: " << diffCorrYCoNum << endl << endl;

Info << "Non-restricting stability numbers:" << endl << endl;

Info<< "Vapor Diffusion Stability Number mean: " << diffYMeanCoNum
    << " max: " << diffYCoNum << endl;

Info<< "Viscosity Stability Number mean: " << muMeanCoNum
    << " max: " << muCoNum << endl << endl;
