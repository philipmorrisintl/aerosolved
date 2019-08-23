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
    max     1E-2;
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

        activityCoeff
        {
            type            Zhang;
            firtSpecieName  Water;
            C1              -0.3049;
            C2              -0.8551;
        }

        heatOfVaporization
        {
            active          false;
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
        diffusion
        {
            type        HirschfelderCurtiss;
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

