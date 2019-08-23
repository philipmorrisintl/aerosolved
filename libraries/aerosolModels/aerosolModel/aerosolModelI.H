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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::aerosolThermo&
Foam::aerosolModel::thermo()
{
    return thermo_;
}


inline const Foam::aerosolThermo&
Foam::aerosolModel::thermo() const
{
    return thermo_;
}


inline const Foam::fvMesh& Foam::aerosolModel::mesh() const
{
    return mesh_;
}


inline Foam::driftFluxModel& Foam::aerosolModel::drift()
{
    return drift_();
}


inline const Foam::driftFluxModel& Foam::aerosolModel::drift() const
{
    return drift_();
}


inline const Foam::compressibleTurbulenceModel&
Foam::aerosolModel::turbulence() const
{
    if (!turbulencePtr_)
    {
        FatalErrorInFunction
            << "turbulencePtr_ is empty. Please use "
            << "aerosolModel::setTurbulence "
            << "(compressibleTurbulenceModel& )"
            << abort(FatalError);
    }

    return *turbulencePtr_;
}


inline const Foam::volScalarField& Foam::aerosolModel::rho() const
{
    return turbulence().rho();
}


inline const Foam::surfaceScalarField& Foam::aerosolModel::phi() const
{
    return turbulence().alphaRhoPhi();
}

inline const Foam::volVectorField& Foam::aerosolModel::U() const
{
    return turbulence().U();
}


inline void Foam::aerosolModel::setTurbulence
(
    compressibleTurbulenceModel& turbModel
)
{
    turbulencePtr_ = &turbModel;
}


inline const Foam::dictionary& Foam::aerosolModel::coeffs() const
{
    return coeffs_;
}

inline const Foam::IOdictionary&
Foam::aerosolModel::outputProperties() const
{
    if (!outputPropertiesPtr_.valid())
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return outputPropertiesPtr_();
}


inline Foam::IOdictionary&
Foam::aerosolModel::outputProperties()
{
    if (!outputPropertiesPtr_.valid())
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return outputPropertiesPtr_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::aerosolModel::mvPhi() const
{
    return mvPhi_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::aerosolModel::mvPhiDrift() const
{
    return mvPhiDrift_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::aerosolModel::mvPhiInertial() const
{
    return mvPhiInertial_();
}

inline const Foam::fv::convectionScheme<Foam::scalar>&
Foam::aerosolModel::mvPhiBrownian() const
{
    return mvPhiBrownian_();
}

inline const Foam::surfaceScalarField& Foam::aerosolModel::phiDrift() const
{
    return phiDrift_;
}

inline const Foam::volSymmTensorField& Foam::aerosolModel::tauDrift() const
{
    return tauDrift_;
}

inline const Foam::surfaceScalarField& Foam::aerosolModel::phiInertial() const
{
    return phiInertial_;
}

inline const Foam::surfaceScalarField& Foam::aerosolModel::phiBrownian() const
{
    return phiBrownian_;
}

inline const Foam::volScalarField& Foam::aerosolModel::DDisp() const
{
    return DDisp_;
}

inline const Foam::PtrList<Foam::volScalarField>&
Foam::aerosolModel::DCont() const
{
    return DCont_;
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::aerosolModel::phiEff() const
{
    return phiEff_;
}

// ************************************************************************* //
