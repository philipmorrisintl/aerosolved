FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1;

vertices
(
    (-VARINLET -VARLHALF -0.5)
    (-VARBHALF -VARLHALF -0.5)
    (VARBHALF -VARLHALF -0.5)
    (VAROUTLET -VARLHALF -0.5)

    (-VARINLET -VARHHALF -0.5)
    (-VARBHALF -VARHHALF -0.5)
    (VARBHALF -VARHHALF -0.5)
    (VAROUTLET -VARHHALF -0.5)

    (-VARINLET VARHHALF -0.5)
    (-VARBHALF VARHHALF -0.5)
    (VARBHALF VARHHALF -0.5)
    (VAROUTLET VARHHALF -0.5)

    (-VARINLET VARLHALF -0.5)
    (-VARBHALF VARLHALF -0.5)
    (VARBHALF VARLHALF -0.5)
    (VAROUTLET VARLHALF -0.5)

    (-VARINLET -VARLHALF 0.5)
    (-VARBHALF -VARLHALF 0.5)
    (VARBHALF -VARLHALF 0.5)
    (VAROUTLET -VARLHALF 0.5)

    (-VARINLET -VARHHALF 0.5)
    (-VARBHALF -VARHHALF 0.5)
    (VARBHALF -VARHHALF 0.5)
    (VAROUTLET -VARHHALF 0.5)

    (-VARINLET VARHHALF 0.5)
    (-VARBHALF VARHHALF 0.5)
    (VARBHALF VARHHALF 0.5)
    (VAROUTLET VARHHALF 0.5)

    (-VARINLET VARLHALF 0.5)
    (-VARBHALF VARLHALF 0.5)
    (VARBHALF VARLHALF 0.5)
    (VAROUTLET VARLHALF 0.5)
);

blocks
(
    hex (0 1 5 4 16 17 21 20) (VARNXI VARNYV 1) simpleGrading (VARG2 VARG2 1)
    hex (1 2 6 5 17 18 22 21) (VARNXB VARNYV 1) simpleGrading (1 VARG2 1)
    hex (2 3 7 6 18 19 23 22) (VARNXO VARNYV 1) simpleGrading (VARG1 VARG2 1)
    hex (4 5 9 8 20 21 25 24) (VARNXI VARNYH 1) simpleGrading (VARG2 1 1)
    hex (6 7 11 10 22 23 27 26) (VARNXO VARNYH 1) simpleGrading (VARG1 1 1)
    hex (8 9 13 12 24 25 29 28) (VARNXI VARNYV 1) simpleGrading (VARG2 VARG1 1)
    hex (9 10 14 13 25 26 30 29) (VARNXB VARNYV 1) simpleGrading (1 VARG1 1)
    hex (10 11 15 14 26 27 31 30) (VARNXO VARNYV 1) simpleGrading (VARG1 VARG1 1)
);

edges
(
);

boundary
(
    lower
    {
        type patch;
        faces
        (
            (16 0 1 17)
            (17 1 2 18)
            (18 2 3 19)
        );
    }
    upper
    {
        type patch;
        faces
        (
            (12 28 29 13)
            (13 29 30 14)
            (14 30 31 15)
        );
    }
    inlet
    {
        type patch;
        faces
        (
            (28 12 8 24)
            (24 8 4 20)
            (20 4 0 16)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (15 31 27 11)
            (11 27 23 7)
            (7 23 19 3)
        );
    }
    wall
    {
        type wall;
        faces
        (
            (21 5 6 22)
            (10 26 22 6)
            (9 25 26 10)
            (25 9 5 21)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (28 24 25 29)
            (24 20 21 25)
            (20 16 17 21)
            (21 17 18 22)
            (22 18 19 23)
            (26 22 23 27)
            (30 26 27 31)
            (29 25 26 30)
            (12 13 9 8)
            (13 14 10 9)
            (14 15 11 10)
            (10 11 7 6)
            (6 7 3 2)
            (5 6 2 1)
            (4 5 1 0)
            (8 9 5 4)
        );
    }
);

mergePatchPairs
(
);
