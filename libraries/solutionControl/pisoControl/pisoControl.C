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

#include "pisoControl.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pisoControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::pisoControl::read()
{
    solutionControl::read(false);

    // Read solution controls
    const dictionary& pisoDict = dict();

    if (!pisoDict.readIfPresent<scalar>("residualCorrector", residualCorrPISO_))
    {
        FatalErrorIn
        (
            "Foam::pisoControl::read"
        )   << "Could not find residualCorrector in solver dictionary" << nl
            << exit(FatalError);
    }

    if (!pisoDict.readIfPresent<scalar>("relTol", relTol_))
    {
        FatalErrorIn
        (
            "Foam::pisoControl::read"
        )   << "Could not find relTol in solver dictionary" << nl
            << exit(FatalError);
    }

    if (!pisoDict.readIfPresent<label>("nCorrectors", nCorrPISO_))
    {
        FatalErrorIn
        (
            "Foam::pisoControl::read"
        )   << "Could not find nCorrectors in solver dictionary" << nl
            << exit(FatalError);
    }

    if (!pisoDict.readIfPresent<label>("nNonOrthogonalCorrectors", nNonOrthCorr_))
    {
        FatalErrorIn
        (
            "Foam::pisoControl::read"
        )   << "Could not find nNonOrthogonalCorrectors in solver dictionary"
            << nl
            << exit(FatalError);
    }

    updateCoeffsPISO_ = pisoDict.lookupOrDefault<Switch>("updateCoeffs", false);
    writeResiduals_ = pisoDict.lookupOrDefault<Switch>("writeResiduals", false);
    theta_ = pisoDict.lookupOrDefault<scalar>("theta", 0.5);
}

Foam::scalar Foam::pisoControl::getResidual()
{
    return residual_;
}

bool Foam::pisoControl::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if (corr_ == 1)
    {
        return false;
    }

    scalar residual = getResidual();

    if (debug)
    {
        Info << "Current residual = " << residual << endl;
    }

    if (residual < residualCorrPISO_)
    {
        if (debug)
        {
            Info << "Convergence criterion met: residual fell below specified tolerance" << endl;;
        }

        return true;
    }
    else if ( (residual/(initPISOresidual_+SMALL)) < relTol_ )
    {
        if (debug)
        {
            Info << "Convergence criterion met: residual fell below specified relative tolerance" << endl;;
        }

        return true;
    }
    else
    {
        return false;
    }
}

void Foam::pisoControl::writeResidual()
{
    // no checks on first iteration - nothing has been calculated yet
    if (corr_ > 1)
    {
        residualsOF_ << mesh_.time().value() << " " << getResidual() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pisoControl::pisoControl(fvMesh& mesh)
:
    solutionControl(mesh, "PISO"),
    nCorrPISO_(0),
    nNonOrthCorr_(0),
    residualCorrPISO_(0.0),
    relTol_(0.0),
    initPISOresidual_(1.0),
    writeResiduals_(false),
    residualsOF_("log.pisoControl"),
    residual_(1.0),
    theta_(0.5),
    InfoLevel_(2)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pisoControl::~pisoControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pisoControl::loop()
{
    read();

    corr_++;

    if (writeResiduals_)
    {
        writeResidual();
    }

    if (corr_ == 2)
    {
        initPISOresidual_ = getResidual();
    }

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    if (corr_ == nCorrPISO_ + 1)
    {
        Info<< algorithmName_ << ": not converged within "
            << nCorrPISO_ << " iterations" << endl;
    }

    // Converged if the criteria are satisfied and if at least two corrections were done

    if (converged())
    {
        corr_ = 0;
        return false;
    }

    return true;
}

bool Foam::pisoControl::converged()
{
    if (corr_ == (nCorrPISO_ + 1) || (criteriaSatisfied() && corr_ > 2))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::pisoControl::updateResidual()
{
    const dictionary& solverDict = mesh_.solverPerformanceDict();
    const List<solverPerformance> sp(solverDict.lookup("p1"));
    residual_ = sp.last().initialResidual();
}

void Foam::pisoControl::silent()
{
    InfoLevel_ = Info.level;
    Info.level = 0;
}

void Foam::pisoControl::verbal()
{
    Info.level = InfoLevel_;
}

// ************************************************************************* //
