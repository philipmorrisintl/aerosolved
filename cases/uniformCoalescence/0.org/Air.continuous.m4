FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Air.continuous;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARY;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
}
