/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (VARX0  VARY0   VARZ0)
    (VARX0  VARY1   VARZ1)
    (VARX0  VARY2   VARZ2)
    (VARX0  VARY3   VARZ3)
    (VARX0  VARY4   VARZ4)
    (VARX0  VARY5   VARZ5)
    (VARX0  VARY6   VARZ6)
    (VARX0  VARY7   VARZ7)

    (VARX1  VARY0   VARZ0)
    (VARX1  VARY1   VARZ1)
    (VARX1  VARY2   VARZ2)
    (VARX1  VARY3   VARZ3)
    (VARX1  VARY4   VARZ4)
    (VARX1  VARY5   VARZ5)
    (VARX1  VARY6   VARZ6)
    (VARX1  VARY7   VARZ7)
);

blocks
(
    hex (4 12 13 5 7 15 14 6) (VARNX VARNBULK VARNBULK) simpleGrading (1 1 1)
    hex (0 8 9 1 4 12 13 5) (VARNX VARNBULK VARNBL) simpleGrading (1 1 VARGBLi)
    hex (5 13 9 1 6 14 10 2) (VARNX VARNBL VARNBULK) simpleGrading (1 VARGBL 1)
    hex (7 15 14 6 3 11 10 2) (VARNX VARNBULK VARNBL) simpleGrading (1 1 VARGBL)
    hex (0 8 12 4 3 11 15 7) (VARNX VARNBL VARNBULK) simpleGrading (1 VARGBLi 1)
);

edges
(
    arc 0 1 (VARX0 VAREY0 VAREZ0)
    arc 1 2 (VARX0 VAREY1 VAREZ1)
    arc 2 3 (VARX0 VAREY2 VAREZ2)
    arc 3 0 (VARX0 VAREY3 VAREZ3)

    arc 4 5 (VARX0 VAREBY0 VAREBZ0)
    arc 5 6 (VARX0 VAREBY1 VAREBZ1)
    arc 6 7 (VARX0 VAREBY2 VAREBZ2)
    arc 7 4 (VARX0 VAREBY3 VAREBZ3)

    arc 8 9 (VARX1 VAREY0 VAREZ0)
    arc 9 10 (VARX1 VAREY1 VAREZ1)
    arc 10 11 (VARX1 VAREY2 VAREZ2)
    arc 11 8 (VARX1 VAREY3 VAREZ3)

    arc 12 13 (VARX1 VAREBY0 VAREBZ0)
    arc 13 14 (VARX1 VAREBY1 VAREBZ1)
    arc 14 15 (VARX1 VAREBY2 VAREBZ2)
    arc 15 12 (VARX1 VAREBY3 VAREBZ3)
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (4 0 1 5)
            (5 1 2 6)
            (3 7 6 2)
            (0 4 7 3)
            (4 5 6 7)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (13 9 8 12)
            (10 9 13 14)
            (10 14 15 11)
            (15 12 8 11)
            (14 13 12 15)
        );
    }

    walls
    {
        type wall;
        faces
        (
            (1 0 8 9)
            (2 1 9 10)
            (3 2 10 11)
            (0 3 11 8)
        );
    }
);

// ************************************************************************* //
