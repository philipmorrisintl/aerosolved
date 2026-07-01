FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      M;
}

dimensions      [-1 0 0 0 0 0 0];

internalField   uniform VARM0;

boundaryField
{
    walls
    {
        type        VARMWALLBCNAME;
        CMD         VARCMD;
        sigma       VARSIGMA;
        gamma       1;
        value       $internalField;
    }
}
