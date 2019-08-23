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
    volScalarField& rDeltaT = trDeltaT.ref();
    volScalarField& rDriftDeltaT = trDriftDeltaT.ref();
    volScalarField& rOrigDeltaT = trOrigDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(readScalar(pimpleDict.lookup("maxCo")));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    // Maximum change in cell concentration per iteration
    // (relative to reference value)
    scalar alphaY(pimpleDict.lookupOrDefault("alphaY", 1.0));

    Info<< "Time scales min/max:" << endl;

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Flow time scale
    {
        rDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phi))()()
           /((2*maxCo)*mesh.V()*rho())
        );

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Flow        = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;
    }

    rOrigDeltaT.ref()=rDeltaT;

    // Drift time scale
    {
        const surfaceScalarField& phiInertial = aerosol->phiInertial();
        const surfaceScalarField& phiBrownian = aerosol->phiBrownian();

        surfaceScalarField phiDrift(phi+phiInertial+phiBrownian);
        rDriftDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phiDrift))()()
            /((2*maxCo)*mesh.V()*rho())
        );
        // Limit the largest time scale
        rDriftDeltaT.max(1/maxDeltaT);
        Info<< "    Drift       = "
            << 1/gMax(rDriftDeltaT.primitiveField()) << ", "
            << 1/gMin(rDriftDeltaT.primitiveField()) << endl;

        rDeltaT.ref() = max(rDeltaT(), rDriftDeltaT());
    }

    // Condensation rate time scale
    if (alphaY < 1)
    {
        dictionary Yref(pimpleDict.subDict("Yref"));

        volScalarField::Internal rDeltaTY
        (
            IOobject
            (
                "rDeltaTY",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("rDeltaTY", rDeltaT.dimensions(), 0)
        );

        bool foundY = false;

        forAll(Z, j)
        {
            if (j != inertIndex)
            {
                volScalarField& Yj = Y[j];

                if (Yref.found(Yj.member()))
                {
                    foundY = true;
                    scalar Yrefj = readScalar(Yref.lookup(Yj.member()));

                    rDeltaTY.field() = max
                    (
                        mag
                        (
                            aerosol->R(Yj)().source()
                           /((Yrefj*alphaY)*(rho*mesh.V()))
                        ),
                        rDeltaTY
                    );
                }
            }
        }

        if (foundY)
        {
            Info<< "    Composition = "
                << 1/(gMax(rDeltaTY.field()) + VSMALL) << ", "
                << 1/(gMin(rDeltaTY.field()) + VSMALL) << endl;

            rDeltaT.ref() = max(rDeltaT(), rDeltaTY);
        }
        else
        {
            IOWarningIn(args.executable().c_str(), Yref)
                << "Cannot find any active species in Yref " << Yref
                << endl;
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "    Overall     = "
        << 1/gMax(rDeltaT.primitiveField())
        << ", " << 1/gMin(rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //
