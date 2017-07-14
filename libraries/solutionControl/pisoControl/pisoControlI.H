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

inline Foam::label Foam::pisoControl::nCorrPISO() const
{
    return nCorrPISO_;
}

inline Foam::label Foam::pisoControl::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline Foam::scalar Foam::pisoControl::residualCorrPISO() const
{
    return residualCorrPISO_;
}

inline Foam::scalar Foam::pisoControl::relTol() const
{
    return relTol_;
}

inline Foam::Switch Foam::pisoControl::updateCoeffsPISO() const
{
    return updateCoeffsPISO_;
}

inline Foam::scalar Foam::pisoControl::theta() const
{
    return theta_;
}

// ************************************************************************* //
