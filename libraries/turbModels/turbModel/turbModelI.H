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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::volScalarField& Foam::turbModel::muTurb() const
{
    return muTurb_;
}

inline const Foam::volScalarField& Foam::turbModel::kTurb() const
{
    return kTurb_;
}

inline const Foam::volScalarField& Foam::turbModel::muEff() const
{
    return muEff_;
}

inline const Foam::volScalarField& Foam::turbModel::kEff() const
{
    return kEff_;
}

inline Foam::fluidThermo& Foam::turbModel::thermo()
{
    return thermo_;
}

inline const Foam::volVectorField& Foam::turbModel::U() const
{
    return U_;
}

inline const Foam::volScalarField& Foam::turbModel::rho() const
{
    return rho_;
}

inline const Foam::volScalarField& Foam::turbModel::delta() const
{
    return delta_;
}

inline const Foam::scalar& Foam::turbModel::deltaCoeff() const
{
    return deltaCoeff_;
}

// ************************************************************************* //
