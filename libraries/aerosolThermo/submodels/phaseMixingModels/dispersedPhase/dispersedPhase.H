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

#ifndef dispersedPhase_H
#define dispersedPhase_H

#include "phaseMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dispersedPhase Declaration
\*---------------------------------------------------------------------------*/

class dispersedPhase
:
    public phaseMixingModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        dispersedPhase(const dispersedPhase&);

        //- Disallow default bitwise assignment
        void operator=(const dispersedPhase&);


public:

    //- Runtime type information
    TypeName("dispersedPhase");


    // Constructors

        //- Construct from thermo model
        dispersedPhase(const aerosolThermo& thermo);


    //- Destructor
    virtual ~dispersedPhase();


    // Member Functions

        //- Compute the mixing
        virtual tmp<volScalarField> mix
        (
            const volScalarField& phiCont,
            const volScalarField& phiDisp
        ) const;

        //- Compute the mixing for a set of cells
        virtual tmp<scalarField> mix
        (
            const scalarField& phiCont,
            const scalarField& phiDisp,
            const labelList& cells
        ) const;

        //- Compute the mixing for a patch
        virtual tmp<scalarField> mix
        (
            const scalarField& phiCont,
            const scalarField& phiDisp,
            const label patchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
