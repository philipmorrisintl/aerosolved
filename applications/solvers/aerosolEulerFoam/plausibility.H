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

\file plausibility.H
\brief Plausibility class definition

This file defines the Plausibility class. It reads certain constraints from the
plausibilityLimits file in the system directory. If a variable does not meet the
plausibility limits, the solution is written to disk and the simulation stopped.
This allows to capture crashes before the simulation explodes. The limits are:

- rhoMin and rhoMax for mixture density \f$\rho\f$
- TMin and TMax for the temperature \f$T\f$
- YMin and YMax for each vapor mass fraction \f$Y_j\f$
- ZMin and ZMax for each liquid mass fraction \f$Y_j\f$

*/

#ifndef PLAUSIBILITY_H
#define PLAUSIBILITY_H

//- Plausibility check of fields, to abort with saving before uncontrolled failure
class Plausibility
{
    private:

        //- Dictionary
        IOdictionary dict_;

        //- Active or not
        const Switch active_;

        //- Thermo object
        const fluidThermo& thermo_;

        //- Limits
        const scalar rhoMin_;
        const scalar rhoMax_;
        const scalar TMin_;
        const scalar TMax_;
        const scalar YMin_;
        const scalar YMax_;
        const scalar ZMin_;
        const scalar ZMax_;

        //- Time object
        Time& runTime_;

    public:

        //- constructor
        Plausibility
        (
            const fluidThermo& thermo,
            Time& runTime,
            const Switch active
        )
        :
            dict_
            (
                IOobject
                (
                    "plausibilityLimits",
                    thermo.mesh().time().system(),
                    thermo.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            ),
            active_(active),
            thermo_(thermo),
            rhoMin_(readScalar(dict_.subDict("limits").lookup("rhoMin"))),
            rhoMax_(readScalar(dict_.subDict("limits").lookup("rhoMax"))),
            TMin_(readScalar(dict_.subDict("limits").lookup("TMin"))),
            TMax_(readScalar(dict_.subDict("limits").lookup("TMax"))),
            YMin_(readScalar(dict_.subDict("limits").lookup("YMin"))),
            YMax_(readScalar(dict_.subDict("limits").lookup("YMax"))),
            ZMin_(readScalar(dict_.subDict("limits").lookup("ZMin"))),
            ZMax_(readScalar(dict_.subDict("limits").lookup("ZMax"))),
            runTime_(runTime)
        {
            if (active_)
            {
                Info << "Plausibility limits check is active. Using:" << endl;
                Info << "    rhoMin = " << rhoMin_ << endl;
                Info << "    rhoMax = " << rhoMax_ << endl;
                Info << "    TMin = " << TMin_ << endl;
                Info << "    TMax = " << TMax_ << endl;
                Info << "    YMin = " << YMin_ << endl;
                Info << "    YMax = " << YMax_ << endl;
                Info << "    ZMin = " << ZMin_ << endl;
                Info << "    ZMax = " << ZMax_ << endl;
            }
            else
            {
                Info << "Plausibility limits check is inactive." << endl;
            }
        }

        //- check for plausibility
        void check()
        {
            bool writeAndEnd = false;

            // density

            const volScalarField& rho = thermo_.rho();

            if (min(rho).value() < rhoMin_)
            {
                Pout << "Plausibility error: min(rho) = "
                     << min(rho).value() << " < " << rhoMin_ << endl;

                writeAndEnd = true;
            }

            if (max(rho).value() > rhoMax_)
            {
                Pout << "Plausibility error: max(rho) = "
                     << max(rho).value() << " > " << rhoMax_ << endl;

                writeAndEnd = true;
            }

            // temperature

            const volScalarField& T = thermo_.T();

            if (min(T).value() < TMin_)
            {
                Pout << "Plausibility error: min(T) = "
                     << min(T).value() << " < " << TMin_ << endl;

                writeAndEnd = true;
            }

            if (max(T).value() > TMax_)
            {
                Pout << "Plausibility error: max(T) = "
                     << max(T).value() << " > " << TMax_ << endl;

                writeAndEnd = true;
            }

            // mass fractions

            forAll(thermo_.species(), j)
            {
                const volScalarField& Y = thermo_.Y()[j];
                const volScalarField& Z = thermo_.Z()[j];

                if (min(Y).value() < YMin_)
                {
                    Pout << "Plausibility error: min(Y) = "
                         << min(Y).value() << " < " << YMin_ << endl;

                    writeAndEnd = true;
                }

                if (max(Y).value() > YMax_)
                {
                    Pout << "Plausibility error: max(Y) = "
                         << max(Y).value() << " > " << YMax_ << endl;

                    writeAndEnd = true;
                }

                if (min(Z).value() < ZMin_)
                {
                    Pout << "Plausibility error: min(Z) = "
                         << min(Z).value() << " < " << ZMin_ << endl;

                    writeAndEnd = true;
                }

                if (max(Z).value() > ZMax_)
                {
                    Pout << "Plausibility error: max(Z) = "
                         << max(Z).value() << " > " << ZMax_ << endl;

                    writeAndEnd = true;
                }
            }

            // if any is out of range, write and abort

            if (writeAndEnd)
            {
                Pout << "Writing fields and aborting." << endl;
                runTime_.writeAndEnd();
            }

        }
};

#endif
