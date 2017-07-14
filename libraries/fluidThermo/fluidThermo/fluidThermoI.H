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

inline const Foam::fvMesh& Foam::fluidThermo::mesh() const
{
    return mesh_;
}

inline Foam::volScalarField& Foam::fluidThermo::T()
{
    return T_;
}

inline const Foam::volScalarField& Foam::fluidThermo::T() const
{
    return T_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::Y()
{
    return Y_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::Y()
const
{
    return Y_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::W()
{
    return W_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::W()
const
{
    return W_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::Z()
{
    return Z_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::Z()
const
{
    return Z_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::X()
{
    return X_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::fluidThermo::X()
const
{
    return X_;
}

inline const Foam::dimensionedScalar& Foam::fluidThermo::R()
const
{
    return R_;
}

inline Foam::volScalarField& Foam::fluidThermo::rho()
{
    return rho_;
}

inline const Foam::volScalarField& Foam::fluidThermo::rho() const
{
    return rho_;
}

inline Foam::volScalarField& Foam::fluidThermo::psi()
{
    return psi_;
}

inline Foam::dimensionedScalar& Foam::fluidThermo::p0()
{
    return p0_;
}

inline const Foam::dimensionedScalar& Foam::fluidThermo::p0() const
{
    return p0_;
}

inline Foam::volScalarField& Foam::fluidThermo::p1()
{
    return p1_;
}

inline const Foam::volScalarField& Foam::fluidThermo::p1() const
{
    return p1_;
}


inline const Foam::volScalarField& Foam::fluidThermo::CvEff() const
{
    return CvEff_;
}

inline const Foam::volScalarField& Foam::fluidThermo::CpEff() const
{
    return CpEff_;
}

inline const Foam::dictionary& Foam::fluidThermo::species() const
{
    return species_;
}

inline const Foam::dictionary& Foam::fluidThermo::speciesPhaseChange() const
{
    return speciesPhaseChange_;
}

inline Foam::label Foam::fluidThermo::nSpecies() const
{
    return species_.size();
}

inline Foam::label Foam::fluidThermo::nSpeciesPhaseChange() const
{
    return speciesPhaseChange_.size();
}

inline Foam::List<Foam::scalar> Foam::fluidThermo::M() const
{
    return M_;
}

inline Foam::List<Foam::scalar> Foam::fluidThermo::Tc() const
{
    return Tc_;
}

inline Foam::List<Foam::Switch> Foam::fluidThermo::phaseChange() const
{
    return phaseChange_;
}

inline Foam::SquareMatrix<Foam::label> Foam::fluidThermo::diffIndex() const
{
    return diffIndex_;
}

// ************************************************************************* //
