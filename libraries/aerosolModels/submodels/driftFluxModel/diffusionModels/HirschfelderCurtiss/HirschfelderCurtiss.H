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

#ifndef HirschfelderCurtiss_H
#define HirschfelderCurtiss_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HirschfelderCurtiss Declaration
\*---------------------------------------------------------------------------*/

class HirschfelderCurtiss
:
    public diffusionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        HirschfelderCurtiss(const HirschfelderCurtiss&);

        //- Disallow default bitwise assignment
        void operator=(const HirschfelderCurtiss&);


public:

    //- Runtime type information
    TypeName("HirschfelderCurtiss");


    // Constructors

        //- Construct from aerosol model
        HirschfelderCurtiss(aerosolModel& aerosol, const dictionary& dict);


    //- Destructor
    virtual ~HirschfelderCurtiss();


    // Member Functions

        //- Compute the diffusivity given a specie index
        virtual tmp<volScalarField> D(const label j);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
