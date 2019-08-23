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
    max     1E-3;
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
        yMin    1E-24;
        yMax    1E-7;
        N       16;
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
        diffusion
        {
            type        none;
        }

        Brownian
        {
            type        StokesEinstein;
        }

        inertial
        {
            type        VARINERTIALMODEL;
            tolerance   1E-6;
            maxIter     3;
            VMax        10.0;
        }
    }
}
