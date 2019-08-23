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

#include "sectionalDistribution.H"
#include "aerosolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sectionalDistribution, 0);
defineRunTimeSelectionTable(sectionalDistribution, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sectionalDistribution::sectionalDistribution
(
    const word& modelType,
    aerosolModel& aerosol,
    const dictionary& dict,
    const dimensionSet& sizeDims
)
:
    PtrList<section>(),
    aerosolSubModelBase(aerosol, dict, typeName, modelType),
    sizeDims_(sizeDims),
    x_(),
    y_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sectionalDistribution::~sectionalDistribution()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

word sectionalDistribution::sectionName(const label i) const
{
    const label d
    (
        Foam::log10
        (
            scalar(this->size())
        )
      + 1
    );

    const Foam::string sectionName(Foam::name(i));

    return std::string(d-sectionName.length(), '0') + sectionName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
