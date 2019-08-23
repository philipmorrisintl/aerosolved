FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties.continuous;
}

thermoType
{
    type            heAerosolRhoThermo;
    mixture         aerosolPhase;
    transport       const;
    thermo          hConst;
    energy          sensibleInternalEnergy;
    equationOfState perfectGas;
    specie          specie;
}

species
{
    Air
    {
        specie
        {
            molWeight       VARW;
        }
        thermodynamics
        {
            Cp              1e3;
            Hf              0;
        }
        transport
        {
            mu              VARMU;
            Pr              1;
        }
        diffusivities
        {
            Air             FullerSchettlerGiddings;
            Water           FullerSchettlerGiddings;
        }
        properties
        {
            Vd              constant 19.7;
        }
    }

    Water
    {
        specie
        {
            molWeight       18.0153;
        }
        thermodynamics
        {
            Cp              1e3;
            Hf              0;
        }
        transport
        {
            mu              1E-5;
            Pr              1;
        }
        diffusivities
        {
            Water           FullerSchettlerGiddings;
        }
        properties
        {
            Vd              constant 20;
            pSat            NSRDS1 (73.649 -7258.2 -7.3037 4.1653e-06 2);
        }
    }
}
