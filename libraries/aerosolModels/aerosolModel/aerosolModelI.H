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

inline Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::M()
{
    return M_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::M() const
{
    return M_;
}

inline Foam::PtrList<Foam::volVectorField>& Foam::aerosolModel::V()
{
    return V_;
}

inline const Foam::PtrList<Foam::volVectorField>& Foam::aerosolModel::V() const
{
    return V_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::S()
{
    return S_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::S() const
{
    return S_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::J()
{
    return J_;
}

inline const Foam::PtrList<Foam::volScalarField>& Foam::aerosolModel::J() const
{
    return J_;
}

inline Foam::PtrList<Foam::surfaceScalarField>& Foam::aerosolModel::phid()
{
    return phid_;
}

inline const Foam::PtrList<Foam::surfaceScalarField>& Foam::aerosolModel::phid() const
{
    return phid_;
}

inline Foam::volScalarField& Foam::aerosolModel::HvapS()
{
    return HvapS_;
}

inline const Foam::volScalarField& Foam::aerosolModel::HvapS() const
{
    return HvapS_;
}

inline const Foam::List<Foam::scalar>& Foam::aerosolModel::y() const
{
    return y_;
}

inline const Foam::List<Foam::scalar>& Foam::aerosolModel::x() const
{
    return x_;
}

inline Foam::label Foam::aerosolModel::P() const
{
    return P_;
}

inline Foam::label Foam::aerosolModel::q() const
{
    return q_;
}

inline Foam::dictionary& Foam::aerosolModel::params()
{
    return params_;
}

inline const Foam::dictionary& Foam::aerosolModel::params() const
{
    return params_;
}

inline const Foam::Switch& Foam::aerosolModel::doMonitors() const
{
    return doMonitors_;
}

inline const Foam::Switch& Foam::aerosolModel::doDrift() const
{
    return doDrift_;
}

inline const Foam::HashPtrTable<Foam::volScalarField>&
Foam::aerosolModel::scalarMonitors() const
{
    return scalarMonitorPtrs_;
}

inline const Foam::HashPtrTable<Foam::volVectorField>&
Foam::aerosolModel::vectorMonitors() const
{
    return vectorMonitorPtrs_;
}

inline bool Foam::aerosolModel::checkUpdateDropDriftVelFields() const
{
    return updateDropDriftVelFields;
}

inline bool Foam::aerosolModel::checkUpdateDropDriftVelField() const
{
    return updateDropDriftVelField;
}

inline bool Foam::aerosolModel::checkGetDropDiffFluxFields() const
{
    return getDropDiffFluxFields;
}

inline bool Foam::aerosolModel::checkGetDropDiffFluxField() const
{
    return getDropDiffFluxField;
}

inline bool Foam::aerosolModel::checkGetCondRateListCell() const
{
    return getCondRateListCell;
}

inline bool Foam::aerosolModel::checkGetCondRateList() const
{
    return getCondRateList;
}

inline bool Foam::aerosolModel::checkGetCondRateField() const
{
    return getCondRateField;
}

inline bool Foam::aerosolModel::checkGetCondRateFields() const
{
    return getCondRateFields;
}

inline bool Foam::aerosolModel::checkGetEtaGammaListCell() const
{
    return getEtaGammaListCell;
}

inline bool Foam::aerosolModel::checkGetEtaGammaList() const
{
    return getEtaGammaList;
}

inline bool Foam::aerosolModel::checkGetCoaRateField() const
{
    return getCoaRateField;
}

inline bool Foam::aerosolModel::checkGetCoaRateCell() const
{
    return getCoaRateCell;
}

inline bool Foam::aerosolModel::checkGetNucFields() const
{
    return getNucFields;
}

inline void Foam::aerosolModel::checkUpdateDropDriftVelFieldsOrExit() const
{
    if (!checkUpdateDropDriftVelFields())
    {
        FatalErrorIn("Foam::aerosolModels::checkUpdateDropDriftVelFieldsOrExit()")
            << "No updateDropDriftVelFields function was set by a drop velocity model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkUpdateDropDriftVelFieldOrExit() const
{
    if (!checkUpdateDropDriftVelField())
    {
        FatalErrorIn("Foam::aerosolModels::checkUpdateDropDriftVelFieldOrExit()")
            << "No updateDropDriftVelField function was set by a drop velocity model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetDropDiffFluxFieldsOrExit() const
{
    if (!checkGetDropDiffFluxFields())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetDropDiffFluxFieldsOrExit()")
            << "No getDropDiffFluxFields function was set by a drop diffusion velocity model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetDropDiffFluxFieldOrExit() const
{
    if (!checkGetDropDiffFluxField())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetDropDiffFluxFieldOrExit()")
            << "No getDropDiffFluxField function was set by a drop velocity model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCondRateListCellOrExit() const
{
    if (!checkGetCondRateListCell())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCondRateListCellOrExit()")
            << "No getCondRateListCell function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCondRateListOrExit() const
{
    if (!checkGetCondRateList())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCondRateListOrExit()")
            << "No getCondRateList function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCondRateFieldOrExit() const
{
    if (!checkGetCondRateField())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCondRateFieldOrExit()")
            << "No getCondRateField function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCondRateFieldsOrExit() const
{
    if (!checkGetCondRateFields())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCondRateFieldsOrExit()")
            << "No getCondRateFields function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetEtaGammaListCellOrExit() const
{
    if (!checkGetEtaGammaListCell())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetEtaGammaListCellOrExit()")
            << "No getEtaGammaListCell function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetEtaGammaListOrExit() const
{
    if (!checkGetEtaGammaList())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetEtaGammaListOrExit()")
            << "No getEtaGammaList function was set by a condensation model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCoaRateFieldOrExit() const
{
    if (!checkGetCoaRateField())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCoaRateFieldOrExit()")
            << "No getCoaRateField function was set by a coalescence model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetCoaRateCellOrExit() const
{
    if (!checkGetCoaRateCell())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCoaRateCellOrExit()")
            << "No getCoaRateCell function was set by a coalescence model." << exit(FatalError);
    }
}

inline void Foam::aerosolModel::checkGetNucFieldsOrExit() const
{
    if (!checkGetNucFields())
    {
        FatalErrorIn("Foam::aerosolModels::checkGetCondRateOrExit()")
            << "No getNucFields function was set by a nucleation model." << exit(FatalError);
    }
}

// ************************************************************************* //
