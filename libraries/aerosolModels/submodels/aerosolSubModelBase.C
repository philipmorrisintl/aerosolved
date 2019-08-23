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

#include "aerosolSubModelBase.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

aerosolSubModelBase::aerosolSubModelBase(aerosolModel& aerosol)
:
    subModelBase(aerosol.outputProperties()),
    aerosol_(aerosol)
{}

aerosolSubModelBase::aerosolSubModelBase
(
    aerosolModel& aerosol,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    subModelBase
    (
        word::null,
        aerosol.outputProperties(),
        dict,
        baseName,
        modelType
    ),
    aerosol_(aerosol)
{}


aerosolSubModelBase::aerosolSubModelBase
(
    const word& modelName,
    aerosolModel& aerosol,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    subModelBase
    (
        modelName,
        aerosol.outputProperties(),
        dict,
        baseName,
        modelType
    ),
    aerosol_(aerosol)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

aerosolSubModelBase::~aerosolSubModelBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool aerosolSubModelBase::writeTime() const
{
    return active() && aerosol_.time().writeTime();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
