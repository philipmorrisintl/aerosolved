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

#include "coalescenceModel.H"
#include "LeeChenLogNormal.H"

#include "makeCoalescenceTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    makeCoalescenceTypes(LeeChenLogNormal, coalescenceModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::LeeChenLogNormal::LeeChenLogNormal
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    coalescenceModel(mesh, aerosol)
{
    if (aerosol.modType() != MOMENTAEROSOLMODEL)
    {
        FatalErrorIn("Foam::coalescenceModels::LeeChenLogNormal::LeeChenLogNormal(const fvMesh& mesh,aerosolModel& aerosol)")
            << "This coalescence model is only valid for a moment method." << exit(FatalError);
    }

    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::LeeChenLogNormal::~LeeChenLogNormal()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::coalescenceModels::LeeChenLogNormal::getCoaRateCell
(
    const Foam::scalar vi,
    const Foam::scalar vj,
    const Foam::label jCell,
    const Foam::scalar Kn
)
{
    const volScalarField& T = thermo_.T();

    const dictionary& species = thermo_.species();
    const label& nSpecies = thermo_.nSpecies();

    const List<Foam::scalar>& Mi = thermo_.M();

    const PtrList<Foam::volScalarField>& Z = thermo_.Z();

    const PtrList<DataEntry<scalar> > dataEntriesListRho_l = thermo_.getProperty("rho", fluidThermo::LIQUID);
    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

    // ------------------------------------------------------------------------
    // Coalescence
    // ------------------------------------------------------------------------

    // Molecular mass
    List<Foam::scalar> m(nSpecies);

    forAll(species, i)
    {
        m[i] = 0.001*Mi[i]/N_A();
    }

    // Liquid density of i-species
    List<Foam::scalar> rho_l(nSpecies);

    forAll(species, i)
    {
        rho_l[i] = dataEntriesListRho_l[i].value(T[jCell]);
    }

    scalar sumZ    = 0.0;
    scalar sumZrho = 0.0;
    forAll(species, i)
    {
        sumZ       += max(Z[i][jCell],0.0);
        sumZrho    += max(Z[i][jCell],0.0)/rho_l[i];
    }

    scalar rho_l_avg = 0.0;
    if (sumZrho > 0.0)
    {
     rho_l_avg = sumZ/sumZrho;
    }

    scalar CMD = pow(vi*6.0/pi(), 1.0/3.0);

    scalar ln2sg = pow(log(W()),2.0);
    scalar K = 2.0*k()*T[jCell]/3.0/muEff[jCell];

    scalar gamma1 = K*(1.0 + exp(ln2sg) + 1.246*Kn*exp(0.5*ln2sg)*(1.0 + exp(2.0*ln2sg) ));

    scalar gamma2 = 0.0;
    if (CMD > SMALL)
    {
        gamma2 = sqrt(3.0*k()*T[jCell]*CMD/rho_l_avg)*pow(1.0 + 1.0/W(),-0.5)*( exp(25.0*ln2sg/8.0) + 2.0*exp(5.0*ln2sg/8.0) + exp(1.0*ln2sg/8.0)  );
    }

    scalar gamma_ =0.0;
    if (gamma1 > 0.0)
    {
      gamma_ = pow(gamma1,-2.0);
    }
    if (gamma2 > 0.0 )
    {
      gamma_ += pow(gamma2,-2.0);
    }

    scalar gamma =0.0;
    if (gamma_ > 0.0)
    {
      gamma = pow(gamma_,-0.5);
    }
    return gamma;
}

bool Foam::coalescenceModels::LeeChenLogNormal::read()
{
    if (coalescenceModel::read())
    {
        // Read coalescence model coefficients

        coeffs_.lookup("W") >> W_;

        thermo_.readProperty("rho", fluidThermo::LIQUID, thermo_.species());

        return true;
    }
    else
    {
        return false;
    }
}

inline const Foam::scalar& Foam::coalescenceModels::LeeChenLogNormal::W() const
{
    return W_;
}

