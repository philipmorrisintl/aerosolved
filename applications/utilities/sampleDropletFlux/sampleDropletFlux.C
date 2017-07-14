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

\file sampleDropletFlux.C
\brief Application to sample the droplet flux at each patch

This application computes the total flux of droplets for each patch and for each
section \f$i\f$. The value that is computed is the numerator of Eq. (7.37) in
\cite thesis, and is given by:

\f[
    \mathrm{total~flux~in~section~}i = \sum_{f\in\mathcal{W}_k} \Phi_{i,f},
\f]

where \f$\mathcal{W}_k\f$ is the set of faces belonging to patch \f$k\f$ and
with

\f[
    \Phi_{i,f} = \left(\phi_f + \phi_{i,f}^{\mathrm{drift}}\right)M_{i,f} +
    \phi_{i,f}^{\mathrm{diff}},
\f]

see Sec. 7.3.3 in \cite thesis.

*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "IOobject.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void updateSizeDistribution
(
    label& P,
    scalarList& x,
    const dictionary& aerosolDict
)
{
    word sizeDistType(aerosolDict.lookup("sizeDistributionType"));

    if (sizeDistType == "none")
    {
        P = 0;
        x.setSize(0);
    }
    else if (sizeDistType == "linear")
    {
        scalar yMin(0.0);
        scalar yMax(0.0);

        aerosolDict.lookup("yMin") >> yMin;
        aerosolDict.lookup("yMax") >> yMax;
        aerosolDict.lookup("P") >> P;

        x.setSize(P);

        forAll(x, i)
        {
            x[i] = yMin + (scalar(i)+0.5)*(yMax-yMin)/P;
        }
    }
    else if (sizeDistType == "logarithmic")
    {
        scalar yMin(0.0);
        scalar yMax(0.0);

        aerosolDict.lookup("yMin") >> yMin;
        aerosolDict.lookup("yMax") >> yMax;
        aerosolDict.lookup("P") >> P;

        x.setSize(P);

        scalar a = Foam::pow(yMax/yMin, 1.0/P);

        x[0] = yMin * Foam::pow(a, 0.5);

        for (label i = 1; i < P; i++)
        {
            x[i] = x[i-1] * a;
        }
    }
    else if (sizeDistType == "geometric")
    {
        scalar yMin(0.0);
        scalar yMax(0.0);
        scalar q(0.0);

        aerosolDict.lookup("yMin") >> yMin;
        aerosolDict.lookup("yMax") >> yMax;
        aerosolDict.lookup("P") >> P;
        aerosolDict.lookup("q") >> q;

        scalarList y(P+1);
        x.setSize(P);

        scalar r = Foam::pow(q, 1.0/(scalar(P)-1.0));

        scalar dy = (yMax - yMin)*(1.0-r) / (1.0-q*r);

        y[0] = yMin;
        y[1] = yMin + dy;

        for (label i = 2; i < P; i++)
        {
            y[i] = y[i-1]*(1+r) - r*y[i-2];
        }

        y[P] = yMax;

        forAll(x, i)
        {
            x[i] = (y[i+1]+r*y[i])/(r+1.0);
        }
    }
    else if (sizeDistType == "list")
    {
        x.setSize(0);

        aerosolDict.lookup("x") >> x;

        P = x.size();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    const dictionary aerosolDict
    (
        IOdictionary
        (
            IOobject
            (
                "aerosolProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict("aerosolModelParameters")
    );

    label P(0);
    scalarList x(0);

    updateSizeDistribution(P, x, aerosolDict);

    Info << "\tx =";

    forAll(x, i)
    {
        Info << " " << x[i];
    }

    Info << endl;

    Info << "\tP = " << P << endl;

    Info << "\tNpatches = " << mesh.boundaryMesh().size() << endl;

    label l(Foam::log10(Foam::scalar(max(P-1,1))) + 1);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "\tTime = " << runTime.timeName() << endl;

        PtrList<surfaceScalarField> phiM(P);

        for(label i = 0; i < P; i++)
        {
            Foam::string is(Foam::name(i));
            Foam::string name("M."+std::string(l-is.length(), '0')+is);

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

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        List< List<scalar> > dropletNumberFlux(bMesh.size(), List<scalar>(P, 0.0));
        Info << "\tPatch names =";

        forAll(bMesh, patchI)
        {
            const polyPatch& p = bMesh[patchI];

            if (p.type() == "patch" || p.type() == "wall")
            {
                forAll(p, localFaceI)
                {
                    for(label i = 0; i < P; i++)
                    {
                        dropletNumberFlux[patchI][i] += phiM[i].boundaryField()[patchI][localFaceI];
                    }
                }

                Info << " " << p.name();
            }
        }
        Info << endl;

        for(label i = 0; i < P; i++)
        {
            Info << "\tSection(" << i << ") =";

            forAll(dropletNumberFlux, patchI)
            {
                Info << " " << dropletNumberFlux[patchI][i];
            }

            Info << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}
