FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (1 0 0);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           $internalField;
    }

    outletFlow
    {
        type            zeroGradient;
    }

    outletProbe
    {
        type            fixedValue;
        value           uniform (VARU 0 0);
    }

    top
    {
        type            slip;
    }

    wedgeFront
    {
        type            wedge;
    }

    wedgeBack
    {
        type            wedge;
    }

    axis
    {
        type            symmetryPlane;
    }

    wallProbe
    {
        type            noSlip;
    }
}
