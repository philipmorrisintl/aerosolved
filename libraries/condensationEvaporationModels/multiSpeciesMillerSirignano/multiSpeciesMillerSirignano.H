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

/**

\file multiSpeciesMillerSirignano.H
\brief Multi-species condensation model based on Miller and Sirignano.

a multi-species condensation model based on Miller et al. \cite Miller:1998aa
and the book of Sirignano \cite Sirignano:1999aa. Three models for the driving
potential function are available: classical rapid mixing, mass analogy I and
mass analogy II, see \cite Miller:1998aa, and two models for the diffusivity
(species-specific and mixture-based). Also two Sherwood numbers can be chosen:
naturation convenction and forced convection.

The model can be configured using three parameters:

- Sherwood = {natural, forced}
- potential = {classical, massAnalogyI, massAnalogyII}
- diffusion = {speciesSpecific, mixture}

*/

#ifndef multiSpeciesMillerSirignano_H
#define multiSpeciesMillerSirignano_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace condensationEvaporationModels
{

enum SherwoodType
{
    NATURALCONVECTION = 0,
    FORCEDCONVECTION = 1
};

enum potentialType
{
    CLASSICAL = 0,
    MASSANALOGYI = 1,
    MASSANALOGYII = 2
};

enum diffusionType
{
    SPECIESSPECIFIC = 0,
    MIXTURE = 1
};

/*---------------------------------------------------------------------------*\
                Class multiSpeciesMillerSirignano Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesMillerSirignano
:
    public condensationEvaporationModel
{
    // Private data

        static const NamedEnum<SherwoodType, 2> SherwoodTypeNames;
        static const NamedEnum<potentialType, 3> potentialTypeNames;
        static const NamedEnum<diffusionType, 2> diffusionTypeNames;


    // Private Member Functions

        //- Disallow copy construct
        multiSpeciesMillerSirignano(const multiSpeciesMillerSirignano&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesMillerSirignano&);


protected:

    // Protected data

        //- Sherwood model
        SherwoodType Sherwood_;

        //- Driving potential model
        potentialType potential_;

        //- Diffusion model
        diffusionType diffusion_;


public:

    //- Runtime type information
    TypeName("multiSpeciesMillerSirignano");


    // Constructors

        //- Construct from components
        multiSpeciesMillerSirignano(const fvMesh& mesh, aerosolModel& aerosol);


    //- Destructor
    virtual ~multiSpeciesMillerSirignano();


    // Member Functions

        // Evolution

            List< List<scalar> > getCondRateListCell
            (
                const List<scalar>& z,
                const label jCell
            );

            List< List<scalar> > getEtaGammaListCell
            (
                const List<scalar>& z,
                const label jCell
            );

            //- Return droplet diameter for a given list of zetas and cell index
            List<scalar> psiInv
            (
                const List<scalar>& zeta,
                const label jCell
            );

        //  Access


        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace condensationEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
