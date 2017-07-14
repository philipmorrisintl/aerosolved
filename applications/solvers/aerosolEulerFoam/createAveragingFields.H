/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2017 Philip Morris International

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

/**

\file createAveragingFields.H
\brief Create fields required for temporal averaging

If the '-timeAveraging' switch is provided, averaging is enabled. In that case
a number of fields will be created globally in this file.

*/

scalar startAveraging = 0;

meanData meanDataFields;

PtrList<volVectorField> meanVectorFields(0);
PtrList<volTensorField> meanTensorFields(0);
PtrList<volScalarField> meanScalarFields(0);
PtrList<volScalarField> meanN(0);
PtrList<volScalarField> meanY(0);
PtrList<volScalarField> meanZ(0);
PtrList<volScalarField> meanS(0);

if(args.options().found("timeAveraging"))
{
    startAveraging = readScalar(runTime.controlDict().lookup("startAveraging"));

    Info << "Starting averaging at " << startAveraging << endl;

    meanVectorFields.setSize(1);
    meanTensorFields.setSize(1);
    meanN.setSize(M.size());
    meanY.setSize(thermo.nSpecies());
    meanZ.setSize(thermo.nSpecies());

    if (aerosol.doMonitors())
    {
        meanScalarFields.setSize(7);
        meanS.setSize(thermo.nSpecies());
    }
    else
    {
        meanScalarFields.setSize(4);
    }

    meanVectorFields.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "mean.U",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("mean.U", dimVelocity, vector::zero)
        )
    );

    meanScalarFields.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "mean.p1",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("mean.p1", dimPressure, 0.0)
        )
    );

    meanTensorFields.set
    (
        0,
        new volTensorField
        (
            IOobject
            (
                "mean.R",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedTensor("mean.R", dimensionSet(0,2,-2,0,0,0,0), tensor::zero)
        )
    );

    meanScalarFields.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "mean.rho",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("mean.rho", dimDensity, 0.0)
        )
    );

    forAll(M, i)
    {
        word name = M[i].name();

        meanN.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    word("mean." + name),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(word("mean." + name), dimless/dimVolume, 0.0)
            )
        );
    }

    forAll(thermo.species(), j)
    {
        word name = thermo.species().keys()[j];

        meanY.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word("mean." + name + "Y"),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(word("mean." + name + "Y"), dimless, 0.0)
            )
        );

        meanZ.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    word("mean." + name + "Z"),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(word("mean." + name + "Z"), dimless, 0.0)
            )
        );
    }

    meanScalarFields.set
    (
        2,
        new volScalarField
        (
            IOobject
            (
                "mean.dcm",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("mean.dcm", dimLength, 0.0)
        )
    );

    meanScalarFields.set
    (
        3,
        new volScalarField
        (
            IOobject
            (
                "mean.dmm",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("mean.dmm", dimLength, 0.0)
        )
    );

    meanDataFields.U = &meanVectorFields[0];
    meanDataFields.p1 = &meanScalarFields[0];
    meanDataFields.R = &meanTensorFields[0];
    meanDataFields.rho = &meanScalarFields[1];
    meanDataFields.N = &meanN;
    meanDataFields.Y = &meanY;
    meanDataFields.Z = &meanZ;
    meanDataFields.dcm = &meanScalarFields[2];
    meanDataFields.dmm = &meanScalarFields[3];

    if (aerosol.doMonitors())
    {
        meanScalarFields.set
        (
            4,
            new volScalarField
            (
                IOobject
                (
                    "mean.Jcoa",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("mean.Jcoa", dimless/dimVolume/dimTime, 0.0)
            )
        );

        meanScalarFields.set
        (
            5,
            new volScalarField
            (
                IOobject
                (
                    "mean.Jnuc",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("mean.Jnuc", dimless/dimVolume/dimTime, 0.0)
            )
        );

        meanScalarFields.set
        (
            6,
            new volScalarField
            (
                IOobject
                (
                    "mean.znuc",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("mean.znuc", dimMass, 0.0)
            )
        );

        forAll(thermo.species(), j)
        {
            const word name("mean.S." + Foam::name(j));

            meanS.set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(name, dimMass/dimVolume/dimTime, 0.0)
                )
            );
        }

        meanDataFields.S = &meanS;
        meanDataFields.Jcoa = &meanScalarFields[4];
        meanDataFields.Jnuc = &meanScalarFields[5];
        meanDataFields.znuc = &meanScalarFields[6];
    }
}
