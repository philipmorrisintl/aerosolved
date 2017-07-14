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

#include "FullerShettlerGiddings.H"
#include "makeDiffusivityModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    makeDiffusivityModels(FullerShettlerGiddings, diffusivityModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::FullerShettlerGiddings::FullerShettlerGiddings
(
    const word& entryName,
    const dictionary& dict,
    const dictionary species,
    const label a,
    const label b
)
:
    diffusivityModel(entryName, species, a, b),
    Vda_(0.0),
    Vdb_(0.0),
    Ma_(0.0),
    Mb_(0.0)
{
    word aName = species_.keys()[a];
    word bName = species_.keys()[b];

    if (species.found(aName) && species.found(bName))
    {
        dictionary av = species.subDict(aName);
        dictionary bv = species.subDict(bName);

        av.subDict("vaporProperties").lookup("Vd") >> Vda_;
        bv.subDict("vaporProperties").lookup("Vd") >> Vdb_;

        av.lookup("moleWeight") >> Ma_;
        bv.lookup("moleWeight") >> Mb_;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::FullerShettlerGiddings::FullerShettlerGiddings(...)"
        )   << "Could not find one of species in dictionary" << nl
            << exit(FatalError);
    }
}


Foam::diffusivityModels::FullerShettlerGiddings::FullerShettlerGiddings
(
    const FullerShettlerGiddings& fsg
)
:
    diffusivityModel(fsg),
    Vda_(fsg.Vda_),
    Vdb_(fsg.Vdb_),
    Ma_(fsg.Ma_),
    Mb_(fsg.Mb_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModels::FullerShettlerGiddings::~FullerShettlerGiddings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::diffusivityModels::FullerShettlerGiddings::value
(
    const Foam::scalar T,
    const Foam::scalar p
) const
{
    return   1E-7
           * pow(T, 1.75)
           / (p/1.01325E+5)
           * sqrt(1.0/Ma_ + 1.0/Mb_)
           / sqr(pow(Vda_, 1.0/3.0) + pow(Vdb_, 1.0/3.0));
}

// ************************************************************************* //
