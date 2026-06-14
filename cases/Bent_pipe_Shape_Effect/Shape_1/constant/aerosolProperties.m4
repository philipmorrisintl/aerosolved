FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      aerosolProperties;
}

aerosolModel    fixedSectional;

diameter
{
    min     1E-9;
    max     1E-3;
}

sfparam
{
   
   particleShape  solid;
}


twoMomentLogNormalAnalyticalCoeffs
{
    sigma   4.0;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    6E-20;
        yMax    9E-12;
        N       25;
    }

    interpolation
    {
        type    twoMoment;
    }

    rescale     false;

    initFromPatch inlet;
}


submodels
{
    condensation
    {
        type        none;
    }

    nucleation
    {
        type        none;
    }

    coalescence
    {
        type        none;
    }

    driftFluxModel
    {
        continuousDiffusion
        {
            type        none;
        }

        dispersedDiffusion
        {
            type        StokesEinstein;
        }

        dispersedInertialDrift
        {
            type        fullStokes;
            tolerance   1E-6;
            maxIter     3;
            VMax        30.0;
            
        }
    }
}


