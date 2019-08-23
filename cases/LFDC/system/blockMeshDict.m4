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
    ( 0.0   0.0 0.0)
    ( 0.75  0.0 0.0)
    ( 0.0  -0.000218096936577626 0.00499524110792016)
    ( 0.75 -0.000218096936577626 0.00499524110792016)
    ( 0.75  0.000218096936577626 0.00499524110792016)
    ( 0.0   0.000218096936577626 0.00499524110792016)
);

blocks
(
    hex (0 1 1 0 2 3 4 5) (VARNX 1 VARNZ) simpleGrading (1 1 1)
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
        );
    }

    inlet
    {
        type patch;
        faces
        (
            (0 2 5 0)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (1 3 4 1)
        );
    }

    walls
    {
        type wall;
        faces
        (
            (2 3 4 5)
        );
    }

    back
    {
        type wedge;
        faces
        (
            (0 1 4 5)
        );
    }

    front
    {
        type wedge;
        faces
        (
            (0 1 3 2)
        );
    }
);

mergePatchPairs
(
);
