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

\file bentPipeDeform.C
\brief Application to bend a mesh

The mesh is always bent around the point \f$\mathbf{x}=(0,0,0)\f$ and around
\f$x_3\f$ such that all \f$x_3\f$ values of all points will remain unchanged.
The radius along which the gemetry is bent is specified by the first command
line argument to this utility.

This application is used to (as the name of the application suggests) deform a
straight pipe into a bent pipe. Such a task could naturally be established by
the blockMesh tool, however, due to a bug (see
https://bugs.openfoam.org/view.php?id=1396), this is difficult to accurately
achieve. With this utility, a perfectly bent geometry can be created.

*/

#include "argList.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("Bend radius");

#   include "setRootCase.H"

    const scalar r = args.argRead<scalar>(1);
    const scalar pi = 3.141592653589793;
    const scalar L = r*pi/2;

#   include "createTime.H"
#   include "createMesh.H"

    // Get times list

    instantList Times = runTime.times();

    pointField zeroPoints(mesh.points());

    // skip "constant" time

    for (label timeI = 1; timeI < Times.size(); ++timeI)
    {
        runTime.setTime(Times[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        pointField newPoints(zeroPoints);

        forAll(zeroPoints, i)
        {
            const point& p0 = zeroPoints[i];
            point& p = newPoints[i];

            if (p0[0] <= 0.0)
            {
                // before the bend, do nothing
            }
            else if (p0[0] > 0.0 && p0[0] <= L)
            {
                // the bend

                const scalar theta = (1.0 - p0[0]/L)*pi/2.0;
                p[0] = Foam::cos(theta) * p0[1];
                p[1] = Foam::sin(theta) * p0[1];
            }
            else
            {
                p[0] = p0[1];
                p[1] = L - p0[0];
            }
        }

        mesh.polyMesh::movePoints(newPoints);
        mesh.write();

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
