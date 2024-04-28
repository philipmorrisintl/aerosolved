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

#include "sectionalLogarithmic.H"
#include "addToRunTimeSelectionTable.H"
#include "aerosolModel.H"
#include "rhoAerosolPhaseThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sectionalLogarithmic, 0);
addToRunTimeSelectionTable
(
    sectionalDistribution,
    sectionalLogarithmic,
    dictionary
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void sectionalLogarithmic::update()
{
    #if OPENFOAM > 1712
    const scalar yMin(dict().get<scalar>("yMin"));
    const scalar yMax(dict().get<scalar>("yMax"));
    const label N(dict().get<label>("N"));
    #else
    const scalar yMin(dict().lookupType<scalar>("yMin"));
    const scalar yMax(dict().lookupType<scalar>("yMax"));
    const label N(dict().lookupType<label>("N"));
    #endif

    if (N < 1)
    {
        FatalErrorInFunction
            << "Minimum sectional distribution size is 1"
            << abort(FatalError);
    }

    if (yMin > yMax)
    {
        FatalErrorInFunction
            << "yMin should be larger than yMax"
            << abort(FatalError);
    }

    if (yMin <= 0)
    {
        FatalErrorInFunction
            << "yMin and yMax should be larger than zero"
            << abort(FatalError);
    }

    // Set x_ and y_ first, so that this information is already accessible

    setSize(N);

    x_.setSize(N);
    y_.setSize(N+1);

    const scalar a(Foam::pow(yMax/yMin, 1.0/N));

    forAll(*this, i)
    {
        x_[i] = yMin*Foam::pow(a, scalar(i)+0.5);
        y_[i] = yMin*Foam::pow(a, scalar(i));
    }

    y_[N] = yMax;

    // Now create the sections

    forAll(*this, i)
    {
        set
        (
            i,
            new section
            (
                i,
                sectionName(i+1),
                yMin*Foam::pow(a, scalar(i)+0.5),
                yMin*Foam::pow(a, scalar(i)),
                yMin*Foam::pow(a, scalar(i)+1.0),
                this->sizeDimensions(),
                aerosol_.mesh()
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sectionalLogarithmic::sectionalLogarithmic
(
    aerosolModel& aerosol,
    const dictionary& dict,
    const dimensionSet& dimensions
)
:
    sectionalDistribution(type(), aerosol, dict, dimensions)
{
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sectionalLogarithmic::~sectionalLogarithmic()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

label sectionalLogarithmic::search(const scalar& s) const
{
    return label
    (
        (Foam::log(s)-Foam::log(yMin()))
      / (Foam::log(yMax())-Foam::log(yMin()))
      * scalar(size())
    );
}

label sectionalLogarithmic::searchLower(const scalar& s) const
{
    return label
    (
        (Foam::log(s)-Foam::log(xMin()))
      / (Foam::log(xMax())-Foam::log(xMin()))
      * scalar(size()-1)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
