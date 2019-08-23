FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1E-3;

vertices
(
    (VARX1 0 0)
    (VARX2 0 0)
    (VARX3 0 0)
    (VARX4 0 0)

    (VARX1 VARY2 -VARZ2)
    (VARX2 VARY2 -VARZ2)
    (VARX3 VARY2 -VARZ2)
    (VARX4 VARY2 -VARZ2)

    (VARX1 VARY3 -VARZ3)
    (VARX2 VARY3 -VARZ3)
    (VARX3 VARY3 -VARZ3)
    (VARX4 VARY3 -VARZ3)

    (VARX1 VARY4 -VARZ4)
    (VARX2 VARY4 -VARZ4)
    (VARX3 VARY4 -VARZ4)
    (VARX4 VARY4 -VARZ4)

    (VARX2 VARY5 -VARZ5)
    (VARX3 VARY5 -VARZ5)
    (VARX4 VARY5 -VARZ5)

    (VARX2 VARY6 -VARZ6)
    (VARX3 VARY6 -VARZ6)

    (VARX1 VARY2 VARZ2)
    (VARX2 VARY2 VARZ2)
    (VARX3 VARY2 VARZ2)
    (VARX4 VARY2 VARZ2)

    (VARX1 VARY3 VARZ3)
    (VARX2 VARY3 VARZ3)
    (VARX3 VARY3 VARZ3)
    (VARX4 VARY3 VARZ3)

    (VARX1 VARY4 VARZ4)
    (VARX2 VARY4 VARZ4)
    (VARX3 VARY4 VARZ4)
    (VARX4 VARY4 VARZ4)

    (VARX2 VARY5 VARZ5)
    (VARX3 VARY5 VARZ5)
    (VARX4 VARY5 VARZ5)

    (VARX2 VARY6 VARZ6)
    (VARX3 VARY6 VARZ6)
);

blocks
(
    hex (0 1 5 4 0 1 22 21) (VARNX1 VARNY1 1) simpleGrading (1 1 1)
    hex (1 2 6 5 1 2 23 22) (VARNX2 VARNY1 1) simpleGrading (VARGX2 1 1)
    hex (2 3 7 6 2 3 24 23) (VARNX3 VARNY1 1) simpleGrading (VARGX3 1 1)
    hex (5 6 10 9 22 23 27 26) (VARNX2 VARNY2 1) simpleGrading (VARGX2 1 1)
    hex (6 7 11 10 23 24 28 27) (VARNX3 VARNY2 1) simpleGrading (VARGX3 1 1)
    hex (8 9 13 12 25 26 30 29) (VARNX1 VARNY3 1) simpleGrading (1 1 1)
    hex (9 10 14 13 26 27 31 30) (VARNX2 VARNY3 1) simpleGrading (VARGX2 1 1)
    hex (10 11 15 14 27 28 32 31) (VARNX3 VARNY3 1) simpleGrading (VARGX3 1 1)
    hex (13 14 17 16 30 31 34 33) (VARNX2 VARNY4 1) simpleGrading (VARGX2 VARGY4 1)
    hex (14 15 18 17 31 32 35 34) (VARNX3 VARNY4 1) simpleGrading (VARGX3 VARGY4 1)
    hex (16 17 20 19 33 34 37 36) (VARNX2 VARNY5 1) simpleGrading (VARGX2 1 1)
);


edges
(
);

boundary
(
    axis
    {
        type empty;
        faces
        (
            (0 1 1 0)
            (1 2 2 1)
            (2 3 3 2)
        );
    }

    inlet1-vapor
    {
        type patch;
        faces
        (
            (21 4 0 0)
        );
    }

    inlet2-heating
    {
        type patch;
        faces
        (
            (29 12 8 25)
        );
    }

    inlet3-cooling
    {
        type patch;
        faces
        (
            (20 19 36 37)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (18 35 32 15)
            (15 32 28 11)
            (11 28 24 7)
            (3 7 24 3)
        );
    }

    walls
    {
        type wall;
        faces
        (
            (4 21 22 5)
            (26 9 5 22)
            (25 8 9 26)
            (12 29 30 13)
            (33 16 13 30)
            (36 19 16 33)
            (20 37 34 17)
            (17 34 35 18)
        );
    }

    back
    {
        type wedge;
        faces
        (
            (0 1 22 21)
            (1 2 23 22)
            (2 3 24 23)
            (22 23 27 26)
            (23 24 28 27)
            (25 26 30 29)
            (26 27 31 30)
            (27 28 32 31)
            (30 31 34 33)
            (31 32 35 34)
            (33 34 37 36)
        );
    }

    front
    {
        type wedge;
        faces
        (
            (0 4 5 1)
            (1 5 6 2)
            (2 6 7 3)
            (5 9 10 6)
            (6 10 11 7)
            (8 12 13 9)
            (9 13 14 10)
            (10 14 15 11)
            (13 16 17 14)
            (14 17 18 15)
            (16 19 20 17)
        );
    }
);

mergePatchPairs
(
);
