FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     aerosolEulerFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         4;

deltaT          1e-6;

writeControl    adjustableRunTime;

writeInterval   0.01;

purgeWrite      0;

writeFormat     ascii;

writePrecision  20;

writeCompression off;

timeFormat      general;

timePrecision   10;

runTimeModifiable yes;

adjustTimeStep  yes;

maxCo           0.5;

maxAerosolCo    0.001;

maxDeltaT       1E-3;

libs            ("libcustomFunctions.so");

functions
{
    dcm
    {
        type            countMeanDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
    }

    dsm
    {
        type            meanDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
        p               3.0;
        q               2.0;
        result          dsm;
    }

    MMD
    {
        type            medianDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
        p               3.0;
        result          MMD;
    }

    Qdot
    {
        type            Qdot;
        lib             ("libaerosolModels.so");
        writeControl    writeTime;
    }

    Cv
    {
        type            thermoField;
        libs            ("libaerosolThermophysicalModels.so");
        writeControl    writeTime;
        thermoField     Cv;
    }

    dropletFlux
    {
        type            VARFLUXTYPE;
        libs            ("libaerosolModels.so");
        patches         (outlet);
        writeControl    writeTime;
    }

    massFlux
    {
        type            massFlux;
        libs            ("libaerosolModels.so");
        patches         (outlet);
        writeControl    writeTime;
    }
}
