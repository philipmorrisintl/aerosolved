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

\file timeAveraging.H
\brief Functions for time averaging

This file contains the definition of the timeAveraging function. A struct of type
meanData must be provided, containing all averaged fields. Based on the specified
startAveraging scalar, which is read from the controlDict in the solver, the
current solution is added to the mean solution using a weight.

*/

#ifndef TIMEAVERAGING_H
#define TIMEAVERAGING_H

struct meanData
{
    volVectorField *U;
    volScalarField *p1, *rho, *dcm, *dmm, *Jcoa, *Jnuc, *znuc;
    volTensorField *R;
    PtrList<volScalarField> *N, *Y, *Z, *S;
};

void timeAveraging
(
    meanData& d,
    const volVectorField& U,
    const fluidThermo& thermo,
    aerosolModel& aerosol,
    const scalar& startAveraging,
    const Time& runTime
)
{
    scalar time = runTime.value();
    scalar dt   = runTime.deltaT().value();

    if (time > startAveraging )
    {
        const volScalarField& p1 = thermo.p1();
        const volScalarField& rho = thermo.rho();

        const PtrList<volScalarField>& M = aerosol.M();

        const PtrList<volScalarField>& Y = thermo.Y();
        const PtrList<volScalarField>& Z = thermo.Z();

        if ((time - dt) < startAveraging)
        {
            // First time

            *d.U == U;
            *d.p1 == p1;
            *d.rho == rho;

            forAll(aerosol.x(), i)
            {
                (*d.N)[i] == M[i]*rho;
            }

            forAll(thermo.species(), j)
            {
                (*d.Y)[j] == Y[j];
                (*d.Z)[j] == Z[j];
            }

            *d.dcm == aerosol.dcm();
            *d.dmm == aerosol.dmm();
        }
        else
        {
            Info << "Averaging is active" <<endl;

            dimensionedScalar tDiff("tDiff", time-startAveraging );
            dimensionedScalar tDiffDt("tDiffDt", tDiff + dt );

            *d.U ==  ((*d.U)*tDiff + U*dt)/tDiffDt;
            *d.p1 == ((*d.p1)*tDiff + p1*dt)/tDiffDt;
            *d.rho == ((*d.rho)*tDiff + rho*dt)/tDiffDt;

            forAll(aerosol.x(), i)
            {
                (*d.N)[i] == ((*d.N)[i]*tDiff + M[i]*rho*dt)/tDiffDt;
            }

            forAll(thermo.species(), j)
            {
                (*d.Y)[j] == ((*d.Y)[j]*tDiff + Y[j]*dt)/tDiffDt;
                (*d.Z)[j] == ((*d.Z)[j]*tDiff + Z[j]*dt)/tDiffDt;
            }

            volTensorField R("R", (U - (*d.U))*(U - (*d.U)));

            *d.R == ((*d.R)*tDiff + R*dt)/tDiffDt;

            if (aerosol.doMonitors())
            {
                const HashPtrTable<volScalarField>& ms
                    = aerosol.scalarMonitors();

                *d.dcm == ((*d.dcm)*tDiff + (*ms["dcm"])*dt)/tDiffDt;
                *d.dmm == ((*d.dmm)*tDiff + (*ms["dmm"])*dt)/tDiffDt;

                *d.Jcoa == ((*d.Jcoa)*tDiff + (*ms["Jcoa"])*dt)/tDiffDt;
                *d.Jnuc == ((*d.Jnuc)*tDiff + (*ms["Jnuc"])*dt)/tDiffDt;
                *d.znuc == ((*d.znuc)*tDiff + (*ms["znuc"])*dt)/tDiffDt;

                forAll(thermo.species(), j)
                {
                    const word name("S." + Foam::name(j));

                    (*d.S)[j] == ((*d.S)[j]*tDiff + (*ms[name])*dt)/tDiffDt;
                }
            }
            else
            {
                *d.dcm == ((*d.dcm)*tDiff + aerosol.dcm()*dt)/tDiffDt;
                *d.dmm == ((*d.dmm)*tDiff + aerosol.dmm()*dt)/tDiffDt;
            }
        }
    }
}

#endif // TIMEAVERAGING_H
