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

/**

\file heAerosolRhoThermo.H
\brief Extension of the heRhoThermo class

This extension of the heRhoThermo class implements the WMix function on top of
standard facilities.

*/

#ifndef heAerosolRhoThermo_H
#define heAerosolRhoThermo_H

#include "heRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class heAerosolRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heAerosolRhoThermo
:
    public heRhoThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        heAerosolRhoThermo
        (
            const heAerosolRhoThermo<BasicPsiThermo, MixtureType>&
        );


public:

    //- Runtime type information
    TypeName("heAerosolRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        heAerosolRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        heAerosolRhoThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~heAerosolRhoThermo();


    // Member functions

        //- Cell-based molecular weight [kg/kmol]
        virtual tmp<scalarField> WMix(const labelList& cells) const;

        //- Patch-based molecular weight [kg/kmol]
        virtual tmp<scalarField> WMix(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heAerosolRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
