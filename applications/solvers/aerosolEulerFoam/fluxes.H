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

\file fluxes.H
\brief Update the mass and number fluxes

This file updates the phiM, phiY, phiZ and phidRho fields. These correspond to
the defitions used in the thesis (see \cite thesis) as follows:

- phiM[i] = \f$\phi_f M_{i,f} + \phi_{i,f}^{\mathrm{drift}} M_{i,f} - \phi_{i,f}^{\mathrm{diff}}\f$,
  see Eq. (7.12).

- phiY[j] = \f$[\phi_f \tilde{Y}_f - \gamma\phi_f^\ell]\frac{Y_{j,f}}{Y_f}\f$, see Eq. (7.18).

- phiZ[j] = \f$[\phi_f \tilde{Z}_f + \phi_f^\ell]\frac{Z_{j,f}}{Z_f}\f$, see Eq. (7.16).

- phidRho = \f$[1-\gamma]\phi_f^\ell\f$, see Eq. (7.20).

*/

linear<scalar> s(mesh);

// Compute gamma

const tmp<volScalarField> tRhol = thermo.rhoLiquid();
const volScalarField& rhol = tRhol();

const tmp<volScalarField> tRhov = thermo.rhoVapor();
const volScalarField& rhov = tRhov();

surfaceScalarField rhovf(s.interpolate(rhov));
surfaceScalarField rholf(s.interpolate(rhol));

thermo.limitVaporDensity(rhovf);
thermo.limitLiquidDensity(rholf);

surfaceScalarField gamma(rhovf/rholf);

// No drift-compensation flux for vapor at walls

forAll(mesh.boundaryMesh(), patchI)
{
    const polyPatch& p = mesh.boundaryMesh()[patchI];

    if (p.type() == "wall")
    {
        gamma.boundaryField()[patchI] = 0.0;
    }
}

// Brownian diffusivities

if (piso.corr() == 0)
{
    brownianDiffusion.update();
}

// Compute fluxes

switch (aerosol.modType())
{

    case MOMENTAEROSOLMODEL:
    {
        const surfaceScalarField phil(phi + phid[0]);
        const surfaceScalarField phiv(phi - gamma*phid[0]);

        tmp<fv::convectionScheme<scalar> > convM
        (
            fv::convectionScheme<scalar>::New(mesh, phil, mesh.divScheme("div(mvLimiter)"))
        );

        phiM[0] = convM->flux(phil, M[0])
                - s.interpolate(DM[0] * rho) * (fvc::snGrad(M[0]) * mesh.magSf());


        // Vapor and liquid mass fraction fluxes

        phidRho = - phi;

        forAll(species, j)
        {
            tmp<fv::convectionScheme<scalar> > convZj
            (
                fv::convectionScheme<scalar>::New(mesh, phil, mesh.divScheme("div(mvLimiter)"))
            );

            phiZ[j] = convZj->flux(phil, Z[j])
                    - s.interpolate(DM[0] * rho) * (fvc::snGrad(Z[j]) * mesh.magSf());

            tmp<fv::convectionScheme<scalar> > convYj
            (
                fv::convectionScheme<scalar>::New(mesh, phiv, mesh.divScheme("div(mvLimiter)"))
            );

            phiY[j] = convYj->flux(phiv, Y[j]);

            phidRho += phiZ[j] + phiY[j];
        }
    }

    break;


    case SECTIONALAEROSOLMODEL:
    {
        surfaceScalarField Ztotf
        (
            IOobject
            (
                "Ztotf",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Ztotf", dimless, 0.0)
        );

        surfaceScalarField Ytotf
        (
            IOobject
            (
                "Ytotf",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Ytotf", dimless, 0.0)
        );

        surfaceScalarField phiSumZ
        (
            IOobject
            (
                "phiSumZ",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("phiSumZ", phi.dimensions(), 0.0)
        );

        surfaceScalarField phidSumZ
        (
            IOobject
            (
                "phidSumZ",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("phidSumZ", phi.dimensions(), 0.0)
        );

        const dimensionedScalar smallDimless("small", dimless, 1E-90);
        const dimensionedScalar unityDimless("unity", dimless, 1.0);


        // Sectional fluxes

        forAll(M, i)
        {
            const dimensionedScalar zi("m", dimMass, aerosol.x()[i]);

            const surfaceScalarField phii(phi + phid[i]);

            surfaceScalarField phidiffi
            (
                - s.interpolate(DM[i] * rho) * (fvc::snGrad(M[i]) * mesh.magSf())
            );

            aerosol.limitWallFlux(phidiffi);

            tmp<fv::convectionScheme<scalar> > convMi
            (
                fv::convectionScheme<scalar>::New(mesh, phii, mesh.divScheme("div(mvLimiter)"))
            );

            phiM[i] = convMi->flux(phii, M[i]) + phidiffi;

            aerosol.limitWallFlux(phiM[i]);

            const surfaceScalarField Mif(convMi->interpolate(phii, M[i]));

            phidSumZ += (phiM[i] - phi * Mif) * zi;

            phiSumZ += phi * Mif * zi;
        }


        // Liquid mass fraction fluxes

        const surfaceScalarField phil(phiSumZ + phidSumZ);

        forAll(species, j)
        {
            tmp<fv::convectionScheme<scalar> > convZj
            (
                fv::convectionScheme<scalar>::New(mesh, phil, mesh.divScheme("div(mvLimiter)"))
            );

            Ztotf += convZj->interpolate(phil, Z[j]);
        }

        forAll(species, j)
        {
            tmp<fv::convectionScheme<scalar> > convZj
            (
                fv::convectionScheme<scalar>::New(mesh, phil, mesh.divScheme("div(mvLimiter)"))
            );

            phiZ[j] = convZj->flux(phil, Z[j]) / stabilise(Ztotf, smallDimless);
        }


        // Vapor mass fraction fluxes

        const surfaceScalarField phiv(phi - phiSumZ - gamma*phidSumZ);

        forAll(species, j)
        {
            tmp<fv::convectionScheme<scalar> > convYj
            (
                fv::convectionScheme<scalar>::New(mesh, phiv, mesh.divScheme("div(mvLimiter)"))
            );

            Ytotf += convYj->interpolate(phiv, Y[j]);
        }

        forAll(species, j)
        {
            tmp<fv::convectionScheme<scalar> > convYj
            (
                fv::convectionScheme<scalar>::New(mesh, phiv, mesh.divScheme("div(mvLimiter)"))
            );

            phiY[j] = convYj->flux(phiv, Y[j]) / stabilise(Ytotf, smallDimless);
        }


        // Density corrections

        phidRho = phil + phiv - phi;
    }

    break;

}
