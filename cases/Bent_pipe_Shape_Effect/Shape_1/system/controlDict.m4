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

endTime         0.5;

deltaT          1e-4;

writeControl    adjustableRunTime;

writeInterval   0.1;

purgeWrite      0;

writeFormat     ascii;

writePrecision  20;

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

        minMax_p
    {
        type          fieldMinMax;
        libs          ("libfieldFunctionObjects.so");
        writeControl  timeStep;
        fields        (p);
    }

    minMax_U
    {
        type          fieldMinMax;
        libs          ("libfieldFunctionObjects.so");
        writeControl  timeStep;
        fields        (U);
    }

   
    minMax_T
    {
        type          fieldMinMax;
        libs          ("libfieldFunctionObjects.so");
        writeControl  timeStep;
        fields        (T);
    }

    yPlus
    {
        type            yPlus;
        libs          ("libfieldFunctionObjects.so");
        patches         (walls);
        writeFields     yes;
        writeControl    writeTime;
    }

     #includeFunc "writeCellCentres"
     #includeFunc "wallShearStress"
     #includeFunc solverInfo

}

