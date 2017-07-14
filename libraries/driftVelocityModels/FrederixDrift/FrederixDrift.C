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

#include "driftVelocityModel.H"
#include "FrederixDrift.H"
#include "fvm.H"
#include "fvc.H"
#include "ddtScheme.H"

#include "makeDriftVelocityTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{
    makeDriftVelocityTypes(FrederixDrift, driftVelocityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModels::FrederixDrift::FrederixDrift
(
    const fvMesh& mesh,
    aerosolModel& aerosol
)
:
    driftVelocityModel(mesh, aerosol),
    g_("g", dimVelocity/dimTime, vector(0, 0, 0)),
    skip_(aerosol.P(), 0),
    r0_(aerosol.P(), 1.0)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModels::FrederixDrift::~FrederixDrift()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::driftVelocityModels::FrederixDrift::updateDropletVelocity
(
    Foam::volVectorField& V,
    const Foam::volScalarField& d,
    const Foam::volScalarField& D,
    const Foam::volVectorField& G,
    const label i
)
{
    V.correctBoundaryConditions();

    label N = mesh_.C().size();
    reduce(N, sumOp<label>());

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    volScalarField DD(D);

    // Iteratively solve the droplet equation of motion for V

    scalar maxRe(0.0);

    const int lvl = Info.level;

    for (label iter = 0; iter < maxIter_; iter++)
    {
        vectorField VPrevIter(V.internalField());

        surfaceScalarField phi(fvc::interpolate(V) & mesh_.Sf());

        // Schiller-Naumann drag coefficient. Else DD is Stokes drag.

        if (SchillerNaumann_)
        {
            tmp<volScalarField> tRe = Re(d, U, V);

            DD.internalField()
                = D.internalField() * (1.0 + 0.15*pow(tRe().internalField(), 0.687));

            maxRe = max(tRe().internalField());
        }

        // The ddt scheme for V.# can be selected by "ddt(V)" in fvSchemes

        fvVectorMatrix VEqn
        (
            fv::ddtScheme<vector>::New
            (
                mesh_,
                mesh_.ddtScheme("ddt(V)")
            )().fvmDdt(V)
          + fvm::div(phi, V, "div(phi,V)")
          - fvm::Sp(fvc::div(phi), V)
          ==
          - fvm::Sp(DD, V) + DD*U
          + G
        );

        // Solve equation (temporarily disable Info output)

        Info.level = 0;

        const scalar r = VEqn.solve(mesh_.solver("V")).initialResidual();

        Info.level = lvl;

        if (iter == 0)
        {
            r0_[i] = r;
        }

        // Compute relative residual

        if (r < TOL_)
        {
            if (!SchillerNaumann_)
            {
                tmp<volScalarField> tRe = Re(d, U, V);
                maxRe = max(tRe().internalField());
            }

            reduce(maxRe, maxOp<scalar>());

            Info << "FrederixDrift: Solving for " << V.name()
                 << ", max(Re) = " << maxRe << ", Initial residual = " << r0_[i]
                 << ", Final residual = " << r
                 << ", No Iterations " << iter+1 << endl;

            break;
        }

        if (iter == (maxIter_-1))
        {
            if (!SchillerNaumann_)
            {
                tmp<volScalarField> tRe = Re(d, U, V);
                maxRe = max(tRe().internalField());
            }

            reduce(maxRe, maxOp<scalar>());

            Info << "FrederixDrift: Field " << V.name()
                 << " droplet velocity equation not converged, max(Re) = "
                 << maxRe << ", Final residual = "
                 << r << " No Iterations " << iter+1 << endl;
        }
    }
}

void Foam::driftVelocityModels::FrederixDrift::updateDropDriftVelFields()
{
    // Must be a sectional aerosol model

    if (aerosol_.modType() != SECTIONALAEROSOLMODEL)
    {
        FatalErrorIn("Foam::driftVelocityModels::FrederixDrift::getDriftVelFields()")
            << "Must be a sectional model" << exit(FatalError);
    }

    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar k("k", dimEnergy/dimTemperature, 1.3806488E-23);

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");
    const volScalarField& T = thermo_.T();

    tmp<volScalarField> tRhol = thermo_.rhoLiquid();
    volScalarField& rhol = tRhol();

    tmp<volScalarField> tRhov = thermo_.rhoVapor();
    volScalarField& rhov = tRhov();

    thermo().limitLiquidDensity(rhol);
    thermo().limitVaporDensity(rhov);

    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const tmp<volScalarField> mg = thermo_.mVapor();

    // Mean free path

    const volScalarField lambda
    (
        sqrt(8.0 * k * T / (pi * mg)) * (4.0/5.0 * muEff / (p0 + p1))
    );

    forAll(aerosol_.x(), i)
    {
        // Check if this one can be skipped

        if (r0_[i] < TOL_)
        {
            if (skip_[i] < maxSkip_)
            {
                skip_[i]++;
                continue;
            }
            else if (skip_[i] == maxSkip_)
            {
                skip_[i] = 0;
            }
        }

        const dimensionedScalar m("z", dimMass, aerosol_.x()[i]);

        const volScalarField d(pow(m/rhol * 6.0/pi, 1.0/3.0));

        const volScalarField Kn(lambda/d);

        volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

        if (!Cunningham_)
        {
            C == dimensionedScalar("one", dimless, 1.0);
        }

        const volScalarField D(18.0 * muEff / (sqr(d) * rhol * C));

        const volVectorField G((rhol-rhov)/rhol * g_);

        updateDropletVelocity(aerosol_.V()[i], d, D, G, i);
    }
}

void Foam::driftVelocityModels::FrederixDrift::updateDropDriftVelField
(
    const Foam::volScalarField& d
)
{
    // Must be a moment aerosol model

    if (aerosol_.modType() != MOMENTAEROSOLMODEL)
    {
        FatalErrorIn("Foam::driftVelocityModels::FrederixDrift::getDriftVelFields()")
            << "Must be a moment model" << exit(FatalError);
    }

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");
    const volScalarField& T = thermo_.T();

    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar k("k", dimEnergy/dimTemperature, 1.3806488E-23);

    tmp<volScalarField> tRhol = thermo_.rhoLiquid();
    volScalarField& rhol = tRhol();

    tmp<volScalarField> tRhov = thermo_.rhoVapor();
    volScalarField& rhov = tRhov();

    thermo().limitLiquidDensity(rhol);
    thermo().limitVaporDensity(rhov);

    const volScalarField& p1 = thermo_.p1();
    const dimensionedScalar p0 = thermo_.p0();

    const tmp<volScalarField> mg = thermo_.mVapor();

    // Mean free path

    const volScalarField lambda
    (
        sqrt(8.0 * k * T / (pi * mg)) * (4.0/5.0 * muEff / (p0 + p1))
    );

    const dimensionedScalar smalld("d", dimLength, 1E-10);

    const volScalarField Kn(lambda/d);

    volScalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

    if (!Cunningham_)
    {
        C == dimensionedScalar("one", dimless, 1.0);
    }

    const volScalarField D(18.0 * muEff / (sqr(stabilise(d, smalld)) * rhol * C));

    const volVectorField G((rhol-rhov)/rhol * g_);

    updateDropletVelocity(aerosol_.V()[0], d, D, G, 0);
}

Foam::tmp<Foam::volScalarField> Foam::driftVelocityModels::FrederixDrift::Re
(
    const Foam::volScalarField& d,
    const Foam::volVectorField& U,
    const Foam::volVectorField& V
)
{
    tmp<volScalarField> tRe
    (
        new volScalarField
        (
            IOobject
            (
                "Re",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Re", dimless, 0.0)
        )
    );

    tmp<volScalarField> tRhov = thermo_.rhoVapor();
    volScalarField& rhov = tRhov();

    thermo().limitVaporDensity(rhov);

    volScalarField& Re = tRe();

    const volScalarField& muEff = mesh_.lookupObject<volScalarField>("muEff");

    Re == d * rhov * mag(U - V) / muEff;

    return tRe;
}

bool Foam::driftVelocityModels::FrederixDrift::read()
{
    if (driftVelocityModel::read())
    {
        params_.lookup("g") >> g_.value();
        params_.lookup("TOL") >> TOL_;
        params_.lookup("maxIter") >> maxIter_;
        params_.lookup("SchillerNaumann") >> SchillerNaumann_;
        params_.lookup("Cunningham") >> Cunningham_;
        params_.lookup("maxSkip") >> maxSkip_;

        if (maxIter_ < 2)
        {
            FatalErrorIn("Foam::driftVelocityModels::FrederixDrift::read()")
                << "maxIter must be at least 2" << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
