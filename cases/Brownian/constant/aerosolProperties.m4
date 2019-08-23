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
    sigma   3;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    1E-19;
        yMax    1E-8;
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
            type        none;
        }
    }
}
