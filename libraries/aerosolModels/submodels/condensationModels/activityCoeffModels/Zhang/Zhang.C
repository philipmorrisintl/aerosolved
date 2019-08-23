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

#include "Zhang.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Zhang, 0);
addToRunTimeSelectionTable(activityCoeffModel, Zhang, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Zhang::Zhang
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    activityCoeffModel(type(), aerosol, dict),
    firstSpecieName_(dict.lookup("firtSpecieName")),
    C1_(readScalar(dict.lookup("C1"))),
    C2_(readScalar(dict.lookup("C2")))
{
    if (aerosol.thermo().activeSpecies().size() != 2)
    {
        FatalErrorInFunction
            << "This activity coefficient model only works for two "
            << "active species" << nl << exit(FatalError);
    }

    if (!aerosol.thermo().activeSpecies().found(firstSpecieName_))
    {
        FatalErrorInFunction
            << "Could not find specie " << firstSpecieName_
            << " in the list of active species" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Zhang::~Zhang()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalarList Zhang::activity(const scalarList& Z) const
{
    aerosolThermo& thermo = aerosol_.thermo();

    rhoAerosolPhaseThermo& thermoCont = thermo.thermoCont();

    const speciesTable& contSpecies = thermo.contSpecies();

    const basicSpecieMixture& compCont = thermoCont.composition();

    const scalar sumZ(min(sum(Z), 1.0));

    scalarList gamma(2, 1.0);

    if (sumZ > VSMALL)
    {
        scalarList W(contSpecies.size(), 0.0);

        forAll(contSpecies, j)
        {
            W[j] = compCont.W(j);
        }

        const scalarList z(Z/sumZ);
        const scalarList w(z/W/sum(z/W));

        const label jA(thermo.activeSpecies()[firstSpecieName_]);
        const label jB(jA == 0 ? 1 : 0);

        gamma[jA] =
            Foam::exp
            (
                C1_/(1.0 + C1_/C2_*(w[jA]/max(1.0-w[jA],VSMALL)))
            );

        gamma[jB] =
            Foam::exp
            (
                C2_/C1_
              * (
                    2.0*Foam::sqrt(C1_*Foam::log(gamma[jA]))
                  + Foam::log(gamma[jA])
                  + C1_
                )
            );
    }

    return gamma;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
