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

#include "blendedCoalescence.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

coaData blendedCoalescence::blend
(
    const coaData& coa1,
    const coaData& coa2,
    const scalar phi
) const
{
    scalarList w(coa1.w()*sqr(phi));
    scalarList p(coa1.p());
    scalarList q(coa1.q());

    w.append(coa2.w()*sqr(1.0-phi));
    p.append(coa2.p());
    q.append(coa2.q());

    return coaData(w, p, q, true);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(blendedCoalescence, 0);
addToRunTimeSelectionTable(coalescenceModel, blendedCoalescence, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blendedCoalescence::blendedCoalescence
(
    aerosolModel& aerosol,
    const dictionary& dict
)
:
    coalescenceModel(type(), aerosol, dict),
    coaModel1_(),
    coaModel2_()
{
    if
    (
        word(dict.subDict("smallKnudsen").lookup("type")) == modelType()
     || word(dict.subDict("largeKnudsen").lookup("type")) == modelType()
    )
    {
        FatalErrorInFunction
            << "Cannot select " << modelType()
            << " as a blended coalescence model"
            << abort(FatalError);
    }

    coaModel1_ = coalescenceModel::New(aerosol, dict.subDict("smallKnudsen"));
    coaModel2_ = coalescenceModel::New(aerosol, dict.subDict("largeKnudsen"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

blendedCoalescence::~blendedCoalescence()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

coaData blendedCoalescence::rate
(
    const scalar& p,
    const scalar& T,
    const scalar& mu,
    const scalar& rhog,
    const scalar& rhol,
    const scalar& d
) const
{
    const coaData coa1(coaModel1_->rate(p, T, mu, rhog, rhol, d));
    const coaData coa2(coaModel2_->rate(p, T, mu, rhog, rhol, d));

    const scalar f(sum(coa1.w()*pow(d, coa1.p()+coa1.q())));
    const scalar g(sum(coa2.w()*pow(d, coa2.p()+coa2.q())));

    // Solution to the equation: phi^2 * f + (1-phi)^2 * g = 1/(1/f+1/g), i.e.,
    // a weight phi is found to reconstruct the harmonic mean of f and g

    const scalar phi(g/(f+g));

    return coaData(blend(coa1, coa2, phi));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
