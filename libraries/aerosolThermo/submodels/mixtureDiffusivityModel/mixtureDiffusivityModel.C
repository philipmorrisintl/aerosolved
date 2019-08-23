/*---------------------------------------------------------------------------*\
License
    AeroSolved
    Copyright (C) 2019 Philip Morris International

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

#include "mixtureDiffusivityModel.H"
#include "aerosolThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * *//

void mixtureDiffusivityModel::getDiffusivityModels()
{
    diffusivities_.clear();
    indices_.clear();

    const rhoAerosolPhaseThermo& thermoCont = thermo_.thermoCont();
    const speciesTable& contSpecies = thermo_.contSpecies();

    const dictionary& speciesDict = thermoCont.speciesDict();

    const label N(contSpecies.size());

    diffusivities_.setSize(label(0.5*N*(N+1)));
    indices_.setSize(N);

    label i = 0;

    forAll(contSpecies, j)
    {
        forAll(contSpecies, k)
        {
            if (j <= k)
            {
                const word jName = contSpecies[j];
                const word kName = contSpecies[k];

                const dictionary& jDict = speciesDict.subDict(jName);
                const dictionary& kDict = speciesDict.subDict(kName);

                if
                (
                    jDict.found("diffusivities") &&
                    jDict.subDict("diffusivities").found(kName)
                )
                {
                    diffusivities_.set
                    (
                        i,
                        diffusivityModel::New
                        (
                            kName,
                            jDict.subDict("diffusivities"),
                            thermo_,
                            j,
                            k
                        )
                    );
                }
                else if
                (
                    kDict.found("diffusivities") &&
                    kDict.subDict("diffusivities").found(jName)
                )
                {
                    diffusivities_.set
                    (
                        i,
                        diffusivityModel::New
                        (
                            jName,
                            kDict.subDict("diffusivities"),
                            thermo_,
                            j,
                            k
                        )
                    );
                }
                else
                {
                    FatalErrorIn
                    (
                        "Foam::mixtureDiffusivityModel::getDiffusivityModels()"
                    )   << "Cannot find binary diffusivity for "
                        << jName << " and " << kName << nl
                        << exit(FatalError);
                }

                indices_[j][k] = i;
                indices_[k][j] = i;

                i++;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixtureDiffusivityModel::mixtureDiffusivityModel
(
    aerosolThermo& thermo
)
:
    thermo_(thermo),
    diffusivities_(),
    indices_()
{
    getDiffusivityModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mixtureDiffusivityModel::~mixtureDiffusivityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

PtrList<scalarField> mixtureDiffusivityModel::Deff() const
{
    const rhoAerosolPhaseThermo& thermoCont = thermo_.thermoCont();
    const speciesTable& contSpecies = thermo_.contSpecies();

    const basicSpecieMixture& compCont = thermoCont.composition();

    scalarList W(contSpecies.size(), 0.0);

    scalarList sumYoverW(thermo_.mesh().nCells(), 0.0);

    forAll(contSpecies, k)
    {
        const scalarField Yk(max(compCont.Y()[k].field(),0.0));

        W[k] = compCont.W(k);

        sumYoverW = sumYoverW + Yk/W[k];
    }

    PtrList<scalarField> Deff(contSpecies.size());

    forAll(Deff, j)
    {
        Deff.set
        (
            j,
            new scalarField(thermo_.mesh().nCells(), 0.0)
        );
    }

    PtrList<scalarField> x(contSpecies.size());

    forAll(x, j)
    {
        const scalarField Yj(max(compCont.Y()[j].field(),0.0));

        x.set(j, new scalarField(Yj/W[j]/sumYoverW));
    }

    forAll(contSpecies, j)
    {
        for (label k = 0; k <= j; k++)
        {
            const scalarList Djk(this->D(j,k));

            Deff[j] += x[k]*Djk;

            if (j != k)
            {
                Deff[k] += x[j]*Djk;
            }
        }
    }

    return Deff;
}


tmp<scalarField> mixtureDiffusivityModel::Deff(const label& j) const
{
    const rhoAerosolPhaseThermo& thermoCont = thermo_.thermoCont();
    const speciesTable& contSpecies = thermo_.contSpecies();

    const basicSpecieMixture& compCont = thermoCont.composition();

    tmp<scalarField> tDeff(new scalarField(thermo_.mesh().nCells(), 0.0));

    scalarField& Deff = tDeff.ref();

    scalarList W(contSpecies.size(), 0.0);

    scalarList sumYoverW(thermo_.mesh().nCells(), 0.0);

    forAll(contSpecies, k)
    {
        const scalarField Yk(max(compCont.Y()[k].field(),0.0));

        W[k] = compCont.W(k);

        sumYoverW = sumYoverW + Yk/W[k];
    }

    forAll(contSpecies, k)
    {
        const scalarField Yk(max(compCont.Y()[k].field(),0.0));

        const scalarField xk(Yk/W[k]/sumYoverW);

        Deff += xk*this->D(j,k);
    }

    return tDeff;
}

tmp<scalarField> mixtureDiffusivityModel::D
(
    const word& jName,
    const word& kName
) const
{
    const speciesTable& contSpecies = thermo_.contSpecies();

    const label& j = contSpecies[jName];
    const label& k = contSpecies[kName];

    return this->D(j,k);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
