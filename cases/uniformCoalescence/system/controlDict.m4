FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    location        "system";
    object          controlDict;
}

application     aerosolEulerFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         VARENDTIME;

deltaT          VARDELTAT;

writeControl    adjustableRunTime;

writeInterval   VARENDTIME;

purgeWrite      0;

writeFormat     binary;

writePrecision  10;

writeCompression no;

timeFormat      general;

timePrecision   12;

runTimeModifiable true;

adjustTimeStep  no;

maxCo           0.5;

maxAerosolCo    0.001;

libs            ("libcustomFunctions.so");

functions
{
    dcm
    {
        type            countMeanDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
    }
    CMD
    {
        type            medianDiameter;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
        p               0.0;
        result          CMD;
    }
    Kn
    {
        type            KnudsenNumber;
        libs            ("libaerosolModels.so");
        writeControl    writeTime;
        result          Kn;
    }
    probes
    {
        type            probes;
        libs            ("libsampling.so");
        writeControl    timeStep;
        writeInterval   1;
        startTime       2E-6;
        probeLocations
        (
            (0.5E-4 0.5E-4 0.5E-4)
        );
        fields
        (
            "(Air|Water|PG|VG).(continuous|dispersed)"
            M
            "M\.[0-9]+"
            dcm
            CMD
            rho
            T
            Kn
        );
    }
}
