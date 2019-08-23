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

#include "coupledCondensation.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledCondensation, 0);
addToRunTimeSelectionTable(condensationModel, coupledCondensation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledCondensation::coupledCondensation
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    condensationModel(type(), aerosol, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coupledCondensation::~coupledCondensation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

conData coupledCondensation::rate
(
    const scalar& p,
    const scalar& T,
    const scalarList& Y,
    const scalarList& Z,
    const scalarList& pSat,
    const scalarList& D,
    const scalarList& rhoCont
) const
{
    // TODO: fix problem with different sizes in condensation/nucleation
    // if(Y.size()!=Z.size()) {
    //     FatalErrorInFunction
    //         << "Sizes of Y " << Y.size() << " and Z "
    //             << Z.size() << " must be equal"
    //             << endl
    //             << exit(FatalError);
    // }

    const scalar pi = constant::mathematical::pi;

    aerosolThermo& thermo = aerosol_.thermo();

    rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

    const basicSpecieMixture& compCont = thermoCont.composition();

    const speciesTable& activeSpecies = thermo.activeSpecies();

    const scalarList Ya(Y, thermo.activeSpeciesMap());

    const scalar sumY(min(sum(Y), 1.0));
    const scalar sumYa(min(sum(Ya), 1.0));
    const scalar sumZ(min(sum(Z), 1.0));

    conData data(activeSpecies.size());

    // Check if we have an adequate mixture

    if (sumZ > 1E-20 && (sumY-sumYa) > 0.0)
    {
        data.active() = true;

        scalarList W(Y.size(), 0.0);

        forAll(Y, j)
        {
            W[j] = compCont.W(j);
        }

        // Kelvin and Fuchs & Sutugin to unity, for now

        scalar Ke = 1.0;
        scalar beta = 1.0;

        // Activity coefficients

        const scalarList gamma(activity_->activity(Z));

        // Compute fractions w.r.t. the dispersed phase

        const scalarList z(Z/sumZ);
        const scalarList w(z/W/sum(z/W));

        // Compute fractions w.r.t. the continuous phase

        const scalarList y(Y/sumY);
        const scalarList x(y/W/sum(y/W));

        // Compute pressures

        const scalarList pSurf(gamma*Ke*pSat*w);
        const scalarList pVapOverY(p/W/sum(Y/W));
        const scalarList pSurfOverZ(gamma*Ke*pSat/W/sum(Z/W));

        // Compute xi

        // Compute diffusivities

        const scalarList xia(x, thermo.inactiveSpeciesMap());
        const scalarList Dia(D, thermo.inactiveSpeciesMap());
        const scalarList Da(D, thermo.activeSpeciesMap());

        const scalar DiaMean(sum(Dia*xia)/sum(xia));

        const scalarList xi
        (
            DiaMean/Da*Foam::log(max(1.0-sum(pSurf)/p, 0.1))
        );

        // Compute condensation rates

        forAll(activeSpecies, j)
        {
            scalar func(0.0);

            if (mag(xi[j]) < 1E-10)
            {
                func = -1.0 + 0.5*xi[j];
            }
            else
            {
                func = xi[j] / (1.0-Foam::exp(xi[j]));
            }

            const scalar c(-2.0*pi*beta*D[j]*rhoCont[j]*func/p);

            data.source()[j] = c*Foam::exp(xi[j])*pVapOverY[j];
            data.sink()[j] = c*pSurfOverZ[j];
        }
    }

    return data;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
