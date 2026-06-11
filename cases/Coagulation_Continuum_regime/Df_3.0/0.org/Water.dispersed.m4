FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Water.dispersed;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARZ;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
}
