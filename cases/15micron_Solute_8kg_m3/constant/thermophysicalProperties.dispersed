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
            rho         1000;
        }
        thermodynamics
        {
            Cp          4195;
            Hf          3.331E5;
        }
        transport
        {
            mu          3.645e-4;
            Pr          2.289;
        }
        properties
        {
            sigma       VDI6 (0.15488 1.64129 -0.75986 -0.85291 1.14113 647.096);
        }
    }
    NaCl
    {
        specie
        {
            molWeight   58.44;
        }
        equationOfState
        {
            rho         2160;
        }
        thermodynamics
        {
            Cp          850;
            Hf          0;
        }
        transport
        {
            mu          1E-3;
            Pr          1;
        }
        properties
        {
            sigma       VDI6 (0.15488 1.64129 -0.75986 -0.85291 1.14113 647.096);
        }
        
    } 
    
}
