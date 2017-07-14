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

#include "condensationEvaporationModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(condensationEvaporationModel, 0);
    defineRunTimeSelectionTable(condensationEvaporationModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::condensationEvaporationModel::condensationEvaporationModel
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    IOdictionary
    (
        IOobject
        (
            "condensationEvaporationProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    aerosol_(aerosol),
    thermo_(aerosol.thermo()),
    params_(subDict("condensationEvaporationModelParameters")),
    I_(0),
    etaGamma_(0)
{
    // Set the aerosol model function pointers to getCondRate()

    aerosol.getCondRateListCell = boost::bind(&condensationEvaporationModel::getCondRateListCell, this, _1, _2);
    aerosol.getCondRateList = boost::bind(&condensationEvaporationModel::getCondRateList, this, _1);
    aerosol.getCondRateField = boost::bind(&condensationEvaporationModel::getCondRateField, this, _1);
    aerosol.getCondRateFields = boost::bind(&condensationEvaporationModel::getCondRateFields, this, _1);
    aerosol.getEtaGammaListCell = boost::bind(&condensationEvaporationModel::getEtaGammaListCell, this, _1, _2);
    aerosol.getEtaGammaList = boost::bind(&condensationEvaporationModel::getEtaGammaList, this, _1);
    aerosol.psiInv = boost::bind(&condensationEvaporationModel::psiInv, this, _1, _2);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::condensationEvaporationModel::~condensationEvaporationModel()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::condensationEvaporationModel::getCondRateList
(
    const Foam::List<Foam::scalar>& z
)
{
    I_.clear();
    I_.setSize(z.size());

    forAll(z, i)
    {
        I_.set(i, new PtrList<volScalarField>(thermo().nSpeciesPhaseChange()));

        forAll(thermo().speciesPhaseChange(), j)
        {
            I_[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        "I",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("I", dimMass/dimTime, 0.0)
                )
            );
        }
    }

    const tmp<volScalarField> tZtot = thermo().Ztot();
    const volScalarField& Ztot = tZtot();

    forAll(mesh_.C(), jCell)
    {
        // Only if we have droplets

        if (Ztot[jCell] > SMALL)
        {
            List<List<scalar> > Icell = getCondRateListCell(z, jCell);

            forAll(z, i)
            {
                forAll(thermo().speciesPhaseChange(), j)
                {
                    PtrList<volScalarField>& Ii = I_[i];

                    volScalarField& I = Ii[j];

                    I[jCell] = Icell[i][j];
                }
            }
        }
    }

    return I_;
}

Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::condensationEvaporationModel::getCondRateFields
(
    const Foam::PtrList<Foam::volScalarField>& z
)
{
    I_.clear();
    I_.setSize(z.size());

    forAll(z, i)
    {
        I_.set(i, new PtrList<volScalarField>(thermo().nSpeciesPhaseChange()));

        forAll(thermo().speciesPhaseChange(), j)
        {
            I_[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        "I",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("I", dimMass/dimTime, 0.0)
                )
            );
        }
    }

    const tmp<volScalarField> tZtot = thermo().Ztot();
    const volScalarField& Ztot = tZtot();

    forAll(mesh_.C(), jCell)
    {
        // Only if we have droplets

        if (Ztot[jCell] > SMALL)
        {
            List<scalar> zList(z.size());

            forAll(z, i)
            {
                zList[i] = z[i][jCell];
            }

            List<List<scalar> > Icell = getCondRateListCell(zList, jCell);

            forAll(z, i)
            {
                forAll(thermo().speciesPhaseChange(), j)
                {
                    I_[i][j][jCell] = Icell[i][j];
                }
            }
        }
    }

    return I_;
}

Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::condensationEvaporationModel::getCondRateField
(
    Foam::volScalarField& z
)
{
    PtrList<volScalarField> listz(1);

    tmp<volScalarField> tz(z);

    listz.set(0, tz);

    return getCondRateFields(listz);
}

Foam::PtrList<Foam::PtrList<Foam::volScalarField> >&
Foam::condensationEvaporationModel::getEtaGammaList
(
    const Foam::List<Foam::scalar>& z
)
{
    etaGamma_.clear();
    etaGamma_.setSize(z.size());

    label n = thermo().nSpeciesPhaseChange();

    forAll(z, i)
    {
        etaGamma_.set(i, new PtrList<volScalarField>(n+1));

        forAll(thermo().speciesPhaseChange(), j)
        {
            etaGamma_[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        word("Gamma." + Foam::name(i) + "." + Foam::name(j)),
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("Gamma", dimless, 0.0)
                )
            );
        }

        etaGamma_[i].set
        (
            n,
            new volScalarField
            (
                IOobject
                (
                    word("eta." + Foam::name(i)),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("eta", dimless, 0.0)
            )
        );
    }

    const tmp<volScalarField> tZtot = thermo().Ztot();
    const volScalarField& Ztot = tZtot();

    forAll(mesh_.C(), jCell)
    {
        // Only if we have droplets

        if (Ztot[jCell] > SMALL)
        {
            List<List<scalar> > etaGammaCell = getEtaGammaListCell(z, jCell);

            forAll(z, i)
            {
                forAll(thermo().speciesPhaseChange(), j)
                {
                    etaGamma_[i][j][jCell] = etaGammaCell[i][j];
                }

                etaGamma_[i][n][jCell] = etaGammaCell[i][n];
            }
        }
    }

    return etaGamma_;
}

bool Foam::condensationEvaporationModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
