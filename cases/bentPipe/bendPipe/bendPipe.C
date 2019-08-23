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

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    const scalar pi(constant::mathematical::pi);

    const scalar R(0.005);
    const scalar RStar(5.7);

    const scalar RBend(R*RStar);
    const scalar D(R*2.0);

    const scalar L(RBend*pi/2.0);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        pointField zeroPoints(mesh.points());

        pointIOField newPoints
            (
                IOobject
                (
                    "points",
                    runTime.findInstance(polyMesh::meshSubDir,"points"),
                    polyMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        forAll(zeroPoints, i)
        {
            const point& p0 = zeroPoints[i];
            point& p = newPoints[i];

            if (p0[0] <= D)
            {
                // before the bend, do nothing
            }
            else if (p0[0] > D && p0[0] <= (D+L))
            {
                // the bend

                const scalar theta = (1.0 - (p0[0]-D)/L)*pi/2.0;

                p[0] = D+Foam::cos(theta) * (p0[1]+RBend);
                p[1] = -RBend+Foam::sin(theta) * (p0[1]+RBend);
            }
            else
            {
                p[0] = D+RBend+p0[1];
                p[1] = -RBend-(p0[0]-(D+L));
            }
        }

        mesh.polyMesh::movePoints(newPoints);
        newPoints.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
