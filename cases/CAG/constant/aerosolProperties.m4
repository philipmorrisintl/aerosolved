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
    min     1E-9;
    max     1E-4;
}

twoMomentLogNormalAnalyticalCoeffs
{
    sigma   1.2;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    1E-25;
        yMax    1E-10;
        N       16;
    }

    interpolation
    {
        type    twoMoment;
    }

    rescale     true;

    initFromPatch inlet1-vapor;
}

submodels
{
    condensation
    {
        type        coupled;

        activityCoeff
        {
            type    constant;
        }

        heatOfVaporization
        {
            active          true;
            limit           true;
            maxDeltaTemp    1.0;
        }
    }

    nucleation
    {
        type        coupled;
        tolerance   1E-7;
    }

    coalescence
    {
        type        none;
    }

    driftFluxModel
    {
        diffusion
        {
            type        none;
        }

        Brownian
        {
            type        none;
        }

        inertial
        {
            type        none;
        }
    }
}
