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
    max     1E-4;
}

fixedSectionalCoeffs
{
    distribution
    {
        type    logarithmic;
        yMin    2E-15;
        yMax    2E-11;
        N       10;
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
            type        constant;
            D           1E-6;
        }

        inertial
        {
            type            VARINERTIALMODEL;

            maxIter         3;
            tolerance       1E-6;

            VMax            10.0;
        }
    }
}
