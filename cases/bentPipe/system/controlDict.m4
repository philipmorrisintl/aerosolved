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

endTime         0.2;

deltaT          1e-4;

writeControl    adjustableRunTime;

writeInterval   0.01;

purgeWrite      0;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   10;

runTimeModifiable yes;

adjustTimeStep  yes;

maxCo           0.5;

maxDeltaT       1E-2;

libs            ("libcustomFunctions.so");

functions
{
    dcm
    {
        type            countMeanDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
    }

    dropletFlux
    {
        type            VARFLUXTYPE;
        libs            ("libaerosolModels.so");
        patches         (inlet outlet walls);
        writeControl    writeTime;
    }

    massFlux
    {
        type            massFlux;
        libs            ("libaerosolModels.so");
        patches         (inlet outlet walls);
        writeControl    writeTime;
    }
}
