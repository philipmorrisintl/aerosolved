FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties.dispersed;
}

thermoType
{
    type            heAerosolRhoThermo;
    mixture         aerosolPhase;
    transport       const;
    thermo          hConst;
    energy          sensibleInternalEnergy;
    equationOfState rhoConst;
    specie          specie;
}

species
{
    Water
    {
        specie
        {
            molWeight   18.015;
        }
        equationOfState
        {
            rho         VARRHOL;
        }
        thermodynamics
        {
            Cp          4195;
            Hf          0;
        }
        transport
        {
            mu          3.645e-4;
            Pr          2.289;
        }
        properties
        {
        }
    }
}
