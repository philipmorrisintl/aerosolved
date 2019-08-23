FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      VG.dispersed;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARVGZ;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
}
