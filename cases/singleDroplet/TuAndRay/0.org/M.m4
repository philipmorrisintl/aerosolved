FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      VARMNAME;
}

dimensions      [-1 0 0 0 0 0 0];

internalField   uniform VARM;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
}
