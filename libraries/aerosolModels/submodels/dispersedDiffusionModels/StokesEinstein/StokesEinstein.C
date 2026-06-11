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

#include "StokesEinstein.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(StokesEinstein, 0);
addToRunTimeSelectionTable(dispersedDiffusionModel, StokesEinstein, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar StokesEinstein::diffusivity
(
    const scalar& d,
    const scalar& p,
    const scalar& T
) const
{
    const scalar& pi = constant::mathematical::pi;
    const scalar& k = constant::physicoChemical::k.value();
    const scalar& NA = constant::physicoChemical::NA.value();

    const rhoAerosolPhaseThermo& thermoCont = aerosol_.thermo().thermoCont();

    const basicSpecieMixture& compCont = thermoCont.composition();

    const label j(thermoCont.species()[aerosol_.thermo().inertSpecie()]);

    const scalar mu(compCont.mu(j, p, T));

    const scalar W(compCont.W(j));

    const scalar mg(0.001*W/NA);

    const scalar lambda
    (
        Foam::sqrt(8.0*k*T/(pi*mg)) * 4.0/5.0*mu/p
    );

    const scalar Kn(2.0*lambda/max(d,aerosol_.dMin()));

    const scalar C(1.0 + (Kn/2.0)*(2.34+1.05*Foam::exp(-0.39/(Kn/2.0))));

    scalar sf = 1.0;
    
    const scalar monoRad = aerosol_.monoRad();
    const scalar pfFm    = aerosol_.dysfPfFm();
    const scalar expFm   = aerosol_.dysfExpFm();
    const scalar pfTr    = aerosol_.dysfPfTr();
    const scalar expTr   = aerosol_.dysfExpTr();
    const scalar pfCont  = aerosol_.dysfPfCont();
    const scalar expCont = aerosol_.dysfExpCont();

	if (Kn < 0.1)
	{
    		sf = pfCont * pow(0.5 * d / monoRad, expCont); // Continuum regime
	}
	else if (Kn > 10.0)
	{
    		sf = pfFm * pow(0.5 * d / monoRad, expFm);     // Free molecular regime
	}
	else
	{
    		sf = pfTr * pow(0.5 * d / monoRad, expTr);     // Transitional regime
	}
      

    return k*T*C/(3.0*pi*mu*d*sf);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StokesEinstein::StokesEinstein
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    dispersedDiffusionModel(type(), aerosol, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

StokesEinstein::~StokesEinstein()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar StokesEinstein::D(const scalar& d, const label& celli) const
{
    const volScalarField& T = aerosol_.thermo().thermoCont().T();
    const volScalarField& p = aerosol_.thermo().p();

    return diffusivity(d, p.field()[celli], T.field()[celli]);
}

scalarField StokesEinstein::D(const scalarField& d, const label& patchi) const
{
    const scalarField& Tp =
        aerosol_.thermo().thermoCont().T().boundaryField()[patchi];

    const scalarField& pp =
        aerosol_.thermo().p().boundaryField()[patchi];

    scalarField Dp(d.size(), 0.0);

    forAll(d, facei)
    {
        Dp[facei] = diffusivity(d[facei], pp[facei], Tp[facei]);
    }

    return Dp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
