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

\file sampleAspirationWedge.C
\brief Application to sample the aerosol aspiration case

This application samples the aerosol aspiration case available in AeroSolved.
This application is specific to that case and cannot be used for anything else.

*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "aerosolModel.H"
#include "fluidThermo.H"
#include "driftVelocityModel.H"
#include "brownianDiffusionModel.H"
#include "nucleationModel.H"
#include "condensationEvaporationModel.H"
#include "coalescenceModel.H"
#include "viscosityModel.H"
#include "List.H"
#include "upwind.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noParallel();
    argList::addOption("inlet", "scalar", "x position of inlet (default 0.0)");
    argList::addOption("probe", "scalar", "x position of probe inlet (default 5*r");
    argList::addOption("outlet", "scalar", "x position of outlet (default 21*r");
    argList::addOption("r", "scalar", "radius of sampling pipe (default 1E-3)");

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    scalar r = args.optionLookupOrDefault("r", 1E-3);

    List<scalar> ps(3, 0.0);

    ps[0] = args.optionLookupOrDefault("inlet", 0.0);
    ps[1] = args.optionLookupOrDefault("probe", 5*r);
    ps[2] = args.optionLookupOrDefault("outlet", 21*r);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        autoPtr<aerosolModel> tAerosolModel(aerosolModel::New(mesh));
        aerosolModel& aerosol = tAerosolModel();

        fluidThermo& thermo = aerosol.thermo();

        autoPtr<viscosityModel> tViscosity(viscosityModel::New(mesh, thermo));
        autoPtr<nucleationModel> tNucleation(nucleationModel::New(mesh, aerosol));
        autoPtr<condensationEvaporationModel> tCondensationEvaporation(condensationEvaporationModel::New(mesh, aerosol));
        autoPtr<coalescenceModel> tCoalescence(coalescenceModel::New(mesh, aerosol));
        autoPtr<driftVelocityModel> tDriftVelocity(driftVelocityModel::New(mesh, aerosol));

        PtrList<volScalarField>& M = aerosol.M();

        PtrList<surfaceScalarField> phiM(M.size());

        forAll(M, i)
        {
            word name = M[i].name();

            phiM.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        word("phi." + name),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        const surfaceScalarField rhof(fvc::interpolate(thermo.rho()));

        List< List<scalar> > dropletMassFlux(aerosol.P(), List<scalar>(3, 0.0));
        List<scalar> volumeFlux(3, 0.0);

        const vectorField& faceCentres = mesh.faceCentres();

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        List<label> pc(3, 0);

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        const volVectorField& C = mesh.C();

        forAll(ps, k)
        {
            forAll(neighbour, faceI)
            {
                const point& faceCentreI = faceCentres[faceI];

                if
                (
                    faceCentreI[0] >= (ps[k]-SMALL) &&
                    faceCentreI[0] <= (ps[k]+SMALL) &&
                    Foam::sqrt(sqr(faceCentreI[1]) + sqr(faceCentreI[2])) <= r
                )
                {
                    label own = owner[faceI];

                    if (C[own][0] > faceCentreI[0])
                    {
                        forAll(aerosol.x(), i)
                        {
                            dropletMassFlux[i][k] -= phiM[i][faceI] * aerosol.x()[i];
                        }

                        volumeFlux[k] -= phi[faceI] / rhof[faceI];
                    }
                    else
                    {
                        forAll(aerosol.x(), i)
                        {
                            dropletMassFlux[i][k] += phiM[i][faceI] * aerosol.x()[i];
                        }

                        volumeFlux[k] += phi[faceI] / rhof[faceI];
                    }
                }
            }

            forAll(bMesh, patchI)
            {
                const polyPatch& p = bMesh[patchI];

                label faceI = p.start();

                forAll(p, localFaceI)
                {
                    const point& faceCentreI = faceCentres[faceI];

                    if
                    (
                        faceCentreI[0] >= (ps[k]-SMALL) &&
                        faceCentreI[0] <= (ps[k]+SMALL) &&
                        Foam::sqrt(sqr(faceCentreI[1]) + sqr(faceCentreI[2])) <= r
                    )
                    {
                        forAll(aerosol.x(), i)
                        {
                            dropletMassFlux[i][k] += phiM[i].boundaryField()[patchI][localFaceI] * aerosol.x()[i];
                        }

                        volumeFlux[k] += phi.boundaryField()[patchI][localFaceI] / rhof.boundaryField()[patchI][localFaceI];

                        ++pc[k];
                    }

                    ++faceI;
                }
            }
        }

        volumeFlux[0] *= -1;

        forAll(aerosol.x(), i)
        {
            dropletMassFlux[i][0] *= -1;
        }

        Info << "Volume fluxes: inlet = " << volumeFlux[0]
             << ", probe = " << volumeFlux[1]
             << ", outlet = " << volumeFlux[2] << endl;

        forAll(aerosol.x(), i)
        {
            Info << "Droplet fluxes, i = " << i << ": inlet = " << dropletMassFlux[i][0]
                 << ", probe = " << dropletMassFlux[i][1]
                 << ", outlet = " << dropletMassFlux[i][2] << endl;
        }

        scalar R = volumeFlux[0]/volumeFlux[2];

        forAll(aerosol.x(), i)
        {
            scalar A(-1.0);

            if (dropletMassFlux[i][0] != 0)
            {
                A = (dropletMassFlux[i][1]/volumeFlux[1]) / (dropletMassFlux[i][0]/volumeFlux[0]);
            }

            Info << "Section: " << runTime.value() << " " << aerosol.x()[i] << " " << " " << R << " " << A << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}
