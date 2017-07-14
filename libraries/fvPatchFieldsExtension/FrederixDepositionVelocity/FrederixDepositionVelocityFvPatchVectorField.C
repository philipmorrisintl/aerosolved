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

#include "FrederixDepositionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "aerosolModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FrederixDepositionVelocityFvPatchVectorField::
FrederixDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    TOL_(0.0),
    Niter_(0),
    Cunningham_(false),
    g_(0, 0, 0),
    SMALL_(1E-10)
{}


Foam::FrederixDepositionVelocityFvPatchVectorField::
FrederixDepositionVelocityFvPatchVectorField
(
    const FrederixDepositionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    TOL_(ptf.TOL_),
    Niter_(ptf.Niter_),
    Cunningham_(ptf.Cunningham_),
    g_(ptf.g_),
    SMALL_(1E-10)
{}


Foam::FrederixDepositionVelocityFvPatchVectorField::
FrederixDepositionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF),
    TOL_(readScalar(dict.lookup("TOL"))),
    Niter_(readLabel(dict.lookup("Niter"))),
    Cunningham_(dict.lookup("Cunningham")),
    g_(dict.lookup("g")),
    SMALL_(1E-10)
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::FrederixDepositionVelocityFvPatchVectorField::
FrederixDepositionVelocityFvPatchVectorField
(
    const FrederixDepositionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(fcvpvf, iF),
    TOL_(fcvpvf.TOL_),
    Niter_(fcvpvf.Niter_),
    Cunningham_(fcvpvf.Cunningham_),
    g_(fcvpvf.g_),
    SMALL_(1E-10)
{
    // Read rhol, for later use

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    thermo.readProperty("rho", fluidThermo::LIQUID, thermo.speciesPhaseChange());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FrederixDepositionVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    zeroGradientFvPatchVectorField::evaluate();

    // Find the necessary references

    const scalar pi = constant::mathematical::pi;

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    const fvPatchField<scalar>& p1 =
        patch().lookupPatchField<volScalarField, scalar>("p1");

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("muEff");

    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>("U");

    // Get a reference to the thermo and aerosol objects

    fluidThermo& thermo =
        const_cast<fluidThermo&>
        (
            db().lookupObject<fluidThermo>("fluidThermoProperties")
        );

    aerosolModel& aerosol =
        const_cast<aerosolModel&>
        (
            db().lookupObject<aerosolModel>("aerosolProperties")
        );

    if (aerosol.modType() != SECTIONALAEROSOLMODEL)
    {
        FatalErrorIn
        (
            "Foam::FrederixDepositionVelocityFvPatchVectorField::evaluate()"
        )   << "This boundary condition only works for a sectional model." << nl
            << exit(FatalError);
    }

    const PtrList<DataEntry<scalar> > dataEntriesListRho_l =
        thermo.getProperty("rho", fluidThermo::LIQUID);

    // Compute liquid density and concentration

    scalarField ZtotRho(patch().size(), 0.0);
    scalarField Ztot(patch().size(), 0.0);
    scalarField rhoZ(patch().size(), 0.0);

    forAll(thermo.species(), j)
    {
        const word name(thermo.species().keys()[j] + "Z");

        scalarField Z
        (
            max(patch().lookupPatchField<volScalarField, scalar>(name), 0.0)
        );

        forAll(ZtotRho, k)
        {
            ZtotRho[k] += Z[k] / stabilise(dataEntriesListRho_l[j].value(T[k]), SMALL);
        }

        Ztot += Z;

        rhoZ += rho*Z;
    }

    scalarField rhol(Ztot/stabilise(ZtotRho, SMALL));

    thermo.limitLiquidDensity(rhol);

    // Patch properties

    vectorField V(*this);

    const vectorField Vc(patchInternalField());
    const vectorField Uc(U.patchInternalField());
    const vectorField n(patch().nf());
    const scalarField delta(mag(this->patch().delta()));

    // Droplet relaxation time

    const label i(sectionNum());

    const scalarField d(pow(aerosol.x()[i]*6.0/pi/rhol, 1.0/3.0));

    scalarField tau(sqr(d)*rhol/(18.0*mu));

    if (Cunningham_)
    {
        const scalar kB(1.3806488E-23);

        const scalar p0 = thermo.p0().value();

        scalarField m(d.size(), 0.0);
        scalarField sumYm(d.size(), 0.0);
        scalarField Ytot(d.size(), 0.0);

        forAll(thermo.species(), j)
        {
            const word name(thermo.species().keys()[j] + "Y");

            const scalar mj = thermo.M()[j]*1E-03/6.022140857E23;

            scalarField Y
            (
                max(patch().lookupPatchField<volScalarField, scalar>(name), 0.0)
            );

            Ytot += Y;

            sumYm += Y/mj;
        }

        Ytot = max(Ytot, 0.0);

        const scalarField mg(Ytot/sumYm);

        const scalarField lambda
        (
             sqrt(8.0 * kB * T / (pi * mg)) * (4.0/5.0 * mu / (p0 + p1))
        );

        const scalarField Kn(lambda/d);

        const scalarField C(1.0 + Kn * (2.34 + 1.05 * exp(-0.39 * 1./Kn)));

        tau *= C;
    }

    // Scaled cell-centered wall-normal fluid and droplet velocity

    const scalarField u((Uc*tau/delta) & n);
    const scalarField v((Vc*tau/delta) & n);

    // Wall-normal scaled gravitational accelleration (gamma is assumed to be small)

    const scalarField g((g_*sqr(tau)/delta) & n);

    forAll(*this, faceI)
    {
        const scalar& uI = u[faceI];
        const scalar& vI = v[faceI];

        const scalar& gI = g[faceI];

        const scalar& deltaI = delta[faceI];
        const scalar& tauI = tau[faceI];
        const vector& nI = n[faceI];

        const scalar r(d[faceI]/deltaI/2.0);

        if (r > 1.0)
        {
            // Cell-center is already within interception distance

            V[faceI] = max(vI, 0.0) * deltaI / tauI * nI;
            continue;
        }

        complex l1, l2;

        if (uI < 0.25)
        {
            l1 = complex(-0.5*sqrt(1.0-4.0*uI)-0.5, 0.0);
            l2 = complex( 0.5*sqrt(1.0-4.0*uI)-0.5, 0.0);
        }
        else
        {
            l1 = complex(-0.5, -0.5*sqrt(4*uI-1.0));
            l2 = complex(-0.5,  0.5*sqrt(4*uI-1.0));
        }

        scalar dzdt(0.0);
        scalar tstar(0.0);

        // Check the mode

        bool noFlowNoG = false;
        bool noFlowFalling1 = false;
        bool noFlowFalling2 = false;
        bool noFlowRising = false;
        bool oscillatory = false;
        bool nonMonotonic = false;
        bool undamped = false;
        bool damped = false;

        if (this->noFlow(l1,l2))
        {
            noFlowNoG = (mag(gI) <= SMALL_ && vI > SMALL_);
            noFlowFalling1 = (gI > SMALL_ && vI > 0.0);
            noFlowFalling2 = (gI > SMALL_ && vI <= 0.0);
            noFlowRising = (gI < -SMALL_ && vI > SMALL_);
        }
        else
        {
            oscillatory = (mag(l1.Im()) > SMALL_);

            nonMonotonic =
            (
                (l1.Re()*vI + (1.0+gI/uI)*l1.Re()*l2.Re()) < (l2.Re()*vI + (1.0+gI/uI)*l1.Re()*l2.Re())
             && (l2.Re()*vI + (1.0+gI/uI)*l1.Re()*l2.Re()) < -SMALL_
             && mag(l1.Re() - l2.Re()) > SMALL_
            );

            undamped =
            (
                (vI > SMALL_ || gI > SMALL_)
             && l2.Re() > SMALL_
             && vI+l1.Re()*(1.0+gI/uI) > SMALL_
            );

            damped =
            (
                l1.Re() < -SMALL_
             && l2.Re() < - SMALL_
             && gI/uI > -r
            );
        }

        if (noFlowNoG)
        {
            // Zero-U and zero-G inertia-induced motion (analytical solution)

            if (vI-1.0 > -r)
            {
                tstar = -log(-(1.0-r)/vI + 1.0);
                dzdt = dxdt(tstar,l1,l2,vI,gI,uI);
            }
        }
        else if (oscillatory || nonMonotonic || undamped || noFlowFalling1 || noFlowFalling2 || noFlowRising)
        {
            // Overshooting motion

            if (noFlowFalling1)
            {
                // Zero-U, positive v, gravity towards the wall

                tstar = (1-r)/min(vI,gI);
            }
            else if (noFlowFalling2)
            {
                // Zero-U, negative or zero v, gravity towards the wall

                tstar = -2.0 * (-gI+vI-1.0+r)/gI;
            }
            else if (noFlowRising)
            {
                // Zero-U, positive v, gravity away from the wall (may not hit the wall)

                tstar = -log(gI/(gI-vI));
            }
            else if (oscillatory)
            {
                // Oscillatory motion

                const scalar a(l1.Re());
                const scalar b(l1.Im());

                const scalar A((a*vI+(1.0+gI/uI)*a*a+(1.0+gI/uI)*b*b)/(b*vI));

                scalarList ts(4, 0.0);

                ts[0] = 2.0*(atan(A+sqrt(A*A+1.0))-1.0*pi)/b;
                ts[1] = 2.0*(atan(A+sqrt(A*A+1.0))+0.0*pi)/b;
                ts[2] = 2.0*(atan(A-sqrt(A*A+1.0))-1.0*pi)/b;
                ts[3] = 2.0*(atan(A-sqrt(A*A+1.0))+0.0*pi)/b;

                scalar zz(-1.0);

                forAll(ts, j)
                {
                    if (ts[j] > 0.0 && x(ts[j],l1,l2,vI,gI,uI) > zz)
                    {
                        zz = x(ts[j],l1,l2,vI,gI,uI);
                        tstar = ts[j];
                    }
                }
            }
            else if (nonMonotonic)
            {
                // Non-monotonic damped or undamped motion

                tstar = - log
                          (
                              (l1.Re()*vI + (1.0+gI/uI)*l1.Re()*l2.Re())
                            / (l2.Re()*vI + (1.0+gI/uI)*l1.Re()*l2.Re())
                          ) / (l1.Re() - l2.Re());
            }
            else if (undamped)
            {
                // Undamped motion

                if (gI/uI > -r)
                {
                    tstar = - log((vI+l2.Re()*(1.0+gI/uI))/(vI+l1.Re()*(1.0+gI/uI)))/(l1.Re()-l2.Re());
                }
                else
                {
                    tstar = 10.0/l2.Re();
                }
            }

            if (tstar > 0.0 && x(tstar,l1,l2,vI,gI,uI) > -r)
            {
                // Droplet hits the wall

                scalar a(0.0);
                scalar b(tstar);

                // Robust bisection

                scalar t(0.0);

                for (label j = 0; j <= Niter_; j++)
                {
                    t = (a+b)/2.0;

                    const scalar z(x(t,l1,l2,vI,gI,uI));

                    if (mag(z+r)/r < TOL_)
                    {
                        break;
                    }

                    if (j == Niter_)
                    {
                        Pout << "Bisection method to find velocity correction did not converge at face "
                             << faceI << ", section = " << i << ", patch "
                             << patch().name() << ". Residual = " << mag(z+r)/r << endl;

                        break;
                    }

                    if (z < -r)
                    {
                        a = t;
                    }
                    else
                    {
                        b = t;
                    }
                }

                dzdt = dxdt(t,l1,l2,vI,gI,uI);
            }
            else
            {
                dzdt = 0.0;
            }
        }
        else if (damped)
        {
            // Non-overshooting damped motion

            tstar = -2.0/l2.Re();

            scalar t = tstar;

            const scalar k(gI/uI);

            for (label j = 0; j <= Niter_; j++)
            {
                scalar z = x(t,l1,l2,vI,gI,uI);
                dzdt = dxdt(t,l1,l2,vI,gI,uI);

                if (mag(z+r)/r < TOL_)
                {
                    break;
                }

                if (j == Niter_)
                {
                    Pout << "Newton method to find velocity correction did not converge at face "
                         << faceI << ", section = " << i << ", patch "
                         << patch().name() << ". Residual = " << mag(z+r)/r << endl;

                    if (!dzdt) dzdt = 0.0;

                    break;
                }

                if ((k-z) < 0.0)
                {
                    t = t/2.0;
                    continue;
                }

                t = t + (k-z)/dzdt * log((k-z)/(k+r));

                if (t < 0.0)
                {
                    // Reduce the initial guess

                    tstar = tstar/2.0;
                    t = tstar;
                }
            }
        }
        else
        {
            dzdt = 0.0;
        }

        V[faceI] = dzdt * deltaI / tauI * nI;
    }

    operator==(V);
}

inline Foam::scalar Foam::FrederixDepositionVelocityFvPatchVectorField::x
(
    const scalar& t,
    const complex& l1,
    const complex& l2,
    const scalar& v,
    const scalar& g,
    const scalar& u
) const
{
    if (!this->noFlow(l1,l2))
    {
        return ((complex(v,0.0)+l2*(1+g/u))/(l1-l2) * expCmplx(t*l1) - (complex(v,0.0)+l1*(1+g/u))/(l1-l2) * expCmplx(t*l2)).Re() + g/u;
    }
    else
    {
        return (g-v)*exp(-t) - (g-v+1.0) + g*t;
    }
}

inline Foam::scalar Foam::FrederixDepositionVelocityFvPatchVectorField::dxdt
(
    const scalar& t,
    const complex& l1,
    const complex& l2,
    const scalar& v,
    const scalar& g,
    const scalar& u
) const
{
    if (!this->noFlow(l1,l2))
    {
        return (l1*(complex(v,0.0)+l2*(1+g/u))/(l1-l2) * expCmplx(t*l1) - l2*(complex(v,0.0)+l1*(1+g/u))/(l1-l2) * expCmplx(t*l2)).Re();
    }
    else
    {
        return -(g-v)*exp(-t) + g;
    }
}

inline bool Foam::FrederixDepositionVelocityFvPatchVectorField::noFlow
(
    const complex& l1,
    const complex& l2
) const
{
    return (mag(l2.Re()) <= SMALL_ && -1-SMALL_ <= l1.Re() && l1.Re() <= -1.0+SMALL_);
}


inline Foam::complex Foam::FrederixDepositionVelocityFvPatchVectorField::expCmplx
(
    const complex& z
) const
{
    return complex
    (
        exp(z.Re())*cos(z.Im()),
        exp(z.Re())*sin(z.Im())
    );
}

inline Foam::word Foam::FrederixDepositionVelocityFvPatchVectorField::name()
{
    return this->dimensionedInternalField().name();
}

inline Foam::label Foam::FrederixDepositionVelocityFvPatchVectorField::sectionNum()
{
    const word name(this->name());

    checkName(name);

    return readLabel(IStringStream(name.substr(2, name.find("_")))());
}

inline bool Foam::FrederixDepositionVelocityFvPatchVectorField::checkName
(
    const Foam::word name
)
{
    regExp r("V\\.([0-9]+)(_[0]+)?");

    if (!r.match(name))
    {
        FatalErrorIn("Foam::sectionalLogNormalFvPatchScalarField::checkName()")
            << "This boundary conditions doens't work on a field named "
            << name
            << exit(FatalError);
    }

    return true;
}

void Foam::FrederixDepositionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("TOL") << TOL_ << token::END_STATEMENT << nl;
    os.writeKeyword("Niter") << Niter_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cunningham") << Cunningham_ << token::END_STATEMENT << nl;
    os.writeKeyword("g") << g_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        FrederixDepositionVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
