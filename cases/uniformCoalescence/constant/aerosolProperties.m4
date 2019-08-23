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
    min     1E-15;
    max     1;
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
        N       VARNSEC;
    }

    interpolation
    {
        type    twoMoment;
    }

    initFromPatch walls;
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
        type        blended;

        smallKnudsen
        {
            type    gasSlip;
            A       1.591;
        }

        largeKnudsen
        {
            type    freeMolecule;
            b       0.70711;
        }
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
