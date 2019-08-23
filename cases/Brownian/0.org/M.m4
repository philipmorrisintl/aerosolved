FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      M;
}

dimensions      [-1 0 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type        VARMBC;
        CMD         constant 1E-8;
        sigma       3;
        value       $internalField;
    }

    outlet
    {
        type        inletOutlet;
        inletValue  $internalField;
        value       $internalField;
    }

    walls
    {
        type        zeroGradient;
    }

    empties
    {
        type        empty;
    }
}
