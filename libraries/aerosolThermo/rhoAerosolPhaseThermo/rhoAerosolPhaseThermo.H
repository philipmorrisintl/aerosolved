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

\file rhoAerosolPhaseThermo.H
\brief Main thermo class for either phase

This class extends the standard rhoReactionThermo class, and adds additional
functionality for the computation of the mole weight and mass density. It is
also based on the aerosolPropertyReader class, which can be used to implement
additional thermo-physical properties (e.g., saturation pressure or surface
tension) into the phase-specific thermo object.

*/

#ifndef rhoAerosolPhaseThermo_H
#define rhoAerosolPhaseThermo_H

#include "rhoReactionThermo.H"
#include "aerosolPropertyReader.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rhoAerosolPhaseThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoAerosolPhaseThermo
:
    public speciesTable,
    public rhoReactionThermo,
    public aerosolPropertyReader
{

public:

    //- Runtime type information
    TypeName("rhoAerosolPhaseThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoAerosolPhaseThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        rhoAerosolPhaseThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoAerosolPhaseThermo> New
        (
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~rhoAerosolPhaseThermo();


    //- Member functions

        inline const speciesTable& species() const
        {
            return *this;
        }

        //- Mixture molecular weight [kg/kmol]
        virtual tmp<volScalarField> WMix() const;

        //- Extend rhoThermo with a cell-based mixture molecular weight
        //- [kg/kmol]. The implementation of this is in heAerosolRhoThermo
        virtual tmp<scalarField> WMix(const labelList& cells) const = 0;

        //- Extend rhoThermo with a patch-based mixture molecular weight
        //- [kg/kmol]. The implementation of this is in heAerosolRhoThermo
        virtual tmp<scalarField> WMix(const label patchi) const = 0;

        //- Density [kg/m3]
        virtual tmp<volScalarField> rho() const;

        //- Extend rhoThermo with a Cell-based density [kg/m3]
        virtual tmp<scalarField> rho(const labelList& cells) const;

        //- Patch-based density [kg/m3]
        virtual tmp<scalarField> rho(const label patchi) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
