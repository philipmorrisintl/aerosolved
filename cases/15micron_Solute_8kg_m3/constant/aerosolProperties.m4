FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      aerosolProperties;
}

aerosolModel    VARMODEL;

diameter
{
    min     1E-12;
    max     1E-2;
}

sfparam
{
   
   particleShape  solid;
}



twoMomentLogNormalAnalyticalCoeffs
{
    sigma   VARSIGMA;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    VARYMIN;
        yMax    VARYMAX;
        N       VARN;
    }

    interpolation
    {
        type    twoMoment;
    }

    rescale     true;
}

submodels
{
    condensation
    {
        type        coupled;
        
        KelvinEffect        true;
        
        FuchsCorrection     true;
        
        DropletTemperatureCorrection true;
        SR    0.55;
        solute       NaCl;
        soluteLimit  0.357;

        heatOfVaporization
        {
            active          true;
        }
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
            type        none;
        }

        dispersedInertialDrift
        {
            type        none;
        }
    }
}

