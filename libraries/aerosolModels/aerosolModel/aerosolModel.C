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

#include "aerosolModel.H"
#include "aerosolModelGitInfo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aerosolModel, 0);
    defineRunTimeSelectionTable(aerosolModel, dictionary);
}

const Foam::word Foam::aerosolModel::aerosolPropertiesName
(
    "aerosolProperties"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::aerosolModel::versionInfo()
{
    Info<< endl;

    Info<< "###############################################################################" << nl
	    << "####                       Welcome to AeroSolved 3.1                       ####" << nl
	    << "###############################################################################" << nl
	    <<  nl
	    << "Cite: " << nl
	    << "Lucci F., Frederix E.M.A., Kuczaj A.K." << nl
	    << "AeroSolved: Computational fluid dynamics modeling of" << nl
	    << "multispecies aerosol flows with sectional and moment methods," << nl
	    << "Journal of Aerosol Science, Volume 159, (2022)" << nl
	    << "doi:10.1016/j.jaerosci.2021.105854" << nl
	    << "Modified Version: 3.1 : Riya Dey, Jayant Krishan, S Anand, Lucci Francesco "<<nl
	    << nl
	    << "###############################################################################" << nl
	    << nl;

    gitInfo();

	Info<< nl
	    << "###############################################################################"
        << nl << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aerosolModel::aerosolModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& aerosolProperties
)
:
    IOdictionary
    (
        IOobject
        (
            aerosolProperties,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    thermo_(mesh),
    turbulencePtr_(),
    mesh_(mesh),
    coeffs_(modelType == "none" ? *this : subDict(modelType + "Coeffs")),
    modelType_(modelType),
    outputPropertiesPtr_(),
    condensation_(),
    nucleation_(),
    coalescence_(),
    drift_(),
    dMin_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("min"))
    ),
    dMax_
    (
        modelType == "none"
      ? 0.0
      : readScalar(subDict("diameter").lookup("max"))
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffs_.lookupOrDefault<scalar>("residualAlpha", 1E-12)
    ),
    
// shape factor parameter

// Particle type: "fractal" or "solid"
   particleShape_
   (
    modelType == "none"
  ? word("none")
  : word(subDict("sfparam").lookup("particleShape"))
   ),


// Set defaults; will overwrite below based on particleShape_
  frDim_(0.0),
  monoRad_(0.0),
  monoRho_(0.0),


// Initialize with defaults (will be set properly below)
	dysfPfFm_(1.0),
	dysfExpFm_(0.0),
	dysfPfTr_(1.0),
	dysfExpTr_(0.0),
	dysfPfCont_(1.0),
	dysfExpCont_(0.0)


// .....................................................       
        
{
      // Validation check here
	if (particleShape_ != "solid" && particleShape_ != "fractal" && particleShape_ != "custom")
	{
    		FatalErrorInFunction
       	 << "Invalid particleShape in sfparam: must be 'solid', 'fractal', or 'custom'. "
       	 << "Found: " << particleShape_ << nl << exit(FatalError);
	}

      const scalar tol = 1e-3;

if (modelType != "none")
{
    
    if (particleShape_ == "solid")
    {
        frDim_       = 3.0;
        monoRad_     = 1.0;  // not used, but set to avoid uninitialized values, should not be non zero
        monoRho_     = 1.0;  // not used, but set to avoid uninitialized values, should not be non zero
        dysfPfFm_    = 1.0;
        dysfExpFm_   = 0.0;
        dysfPfTr_    = 1.0;
        dysfExpTr_   = 0.0;
        dysfPfCont_  = 1.0;
        dysfExpCont_ = 0.0;
    }
    else if (particleShape_ == "fractal")
    
    {       
            frDim_    = readScalar(subDict("sfparam").lookup("frdim"));
	     monoRad_ = readScalar(subDict("sfparam").lookup("monorad"));
             monoRho_ = readScalar(subDict("sfparam").lookup("monorho"));
            
        if (mag(frDim_ - 3.0) < tol)
        {
            dysfPfFm_    = 1.0;
            dysfExpFm_   = 0.0;
            dysfPfTr_    = 1.0;
            dysfExpTr_   = 0.0;
            dysfPfCont_  = 1.0;
            dysfExpCont_ = 0.0;
        }
        else if (mag(frDim_ - 2.49) < tol)
        {
            dysfPfFm_    = 0.93;
            dysfExpFm_   = 0.49;
            dysfPfTr_    = 1.12;
            dysfExpTr_   = 0.28;
            dysfPfCont_  = 0.57;
            dysfExpCont_ = 0.40;
        }
        else if (mag(frDim_ - 2.25) < tol)
        {
            dysfPfFm_    = 0.91;
            dysfExpFm_   = 0.54;
            dysfPfTr_    = 1.24;
            dysfExpTr_   = 0.28;
            dysfPfCont_  = 0.58;
            dysfExpCont_ = 0.43;
        }
        else if (mag(frDim_ - 2.00) < tol)
        {
            dysfPfFm_    = 0.90;
            dysfExpFm_   = 0.56;
            dysfPfTr_    = 1.24;
            dysfExpTr_   = 0.31;
            dysfPfCont_  = 0.76;
            dysfExpCont_ = 0.40;
        }
        else if (mag(frDim_ - 1.80) < tol)
        {
            dysfPfFm_    = 0.90;
            dysfExpFm_   = 0.58;
            dysfPfTr_    = 1.19;
            dysfExpTr_   = 0.36;
            dysfPfCont_  = 0.70;
            dysfExpCont_ = 0.45;
        }
        else
        {
         
         Info<< nl
           << "ERROR: Invalid Df specified." << nl
           << "Df = " << frDim_ << nl
           << "Supported default values for Df in fractal mode are:" << nl
           << "3.0, 2.49, 2.25, 2.00, 1.80" << nl
           << "Rather use Custom mode in aerosolProperties."
           << endl;

         
         FatalErrorIn("aerosolModel::aerosolModel")
        << "Df = " << frDim_ << nl
        << "Df out of default value."
        << nl << "Supported default values for Df in fractal mode are:"
        << nl << "3.0, 2.49, 2.25, 2.00, 1.80"
        << nl <<"Rather use Custom mode in aerosolProperties."
        
        << exit(FatalError);
        }
    }
    
    else if (particleShape_ == "custom")
    {
        frDim_       = readScalar(subDict("sfparam").lookup("frdim"));
        monoRad_     = readScalar(subDict("sfparam").lookup("monorad"));
        monoRho_     = readScalar(subDict("sfparam").lookup("monorho"));
        dysfPfFm_    = readScalar(subDict("sfparam").lookup("dysfPfFm"));
        dysfExpFm_   = readScalar(subDict("sfparam").lookup("dysfExpFm"));
        dysfPfTr_    = readScalar(subDict("sfparam").lookup("dysfPfTr"));
        dysfExpTr_   = readScalar(subDict("sfparam").lookup("dysfExpTr"));
        dysfPfCont_  = readScalar(subDict("sfparam").lookup("dysfPfCont"));
        dysfExpCont_ = readScalar(subDict("sfparam").lookup("dysfExpCont"));
    }

}

	Info << "\n--- Shape Factor Configuration ---" << nl;
	Info << "particleShape     : " << particleShape_ << nl;
	Info << "frDim            : " << frDim_ << nl;
	Info << "monorad          : " << monoRad_ << nl;
	Info << "monorho          : " << monoRho_ << nl;
	Info << "dysfPfFm         : " << dysfPfFm_ << nl;
	Info << "dysfExpFm        : " << dysfExpFm_ << nl;
	Info << "dysfPfTr         : " << dysfPfTr_ << nl;
	Info << "dysfExpTr        : " << dysfExpTr_ << nl;
	Info << "dysfPfCont       : " << dysfPfCont_ << nl;
	Info << "dysfExpCont      : " << dysfExpCont_ << nl;
	Info << "----------------------------------\n" << endl;
    
    
    

    versionInfo();

    read();

    if (!outputPropertiesPtr_.valid())
    {
        const fileName uniformPath(word("uniform")/"aerosolModels");

        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "outputProperties",
                    mesh_.time().timeName(),
                    uniformPath,
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    if (modelType != "none")
    {
        condensation_ =
            condensationModel::New
            (
                *this,
                subDict("submodels").subDict("condensation")
            );

        nucleation_ =
            nucleationModel::New
            (
                *this,
                subDict("submodels").subDict("nucleation")
            );

        coalescence_ =
            coalescenceModel::New
            (
                *this,
                subDict("submodels").subDict("coalescence")
            );
    }

    drift_.reset(new driftFluxModel(*this, subDict("submodels")));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aerosolModel::~aerosolModel()
{
    if (turbulencePtr_)
    {
        turbulencePtr_ = 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::aerosolModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}

Foam::tmp<Foam::scalarField> Foam::aerosolModel::getRDeltaT()
{
    if(mesh_.foundObject<volScalarField>("rDeltaT"))
    {
        return tmp<scalarField>
        (
            mesh_.lookupObject<volScalarField>("rDeltaT")
        );
    }
    else
    {
        const scalarField& rho = this->rho().field();

        if(!rDeltaT_.valid())
        {
            rDeltaT_.reset(new scalarField(rho.size()));
        }

        rDeltaT_() = 1.0/mesh_.time().deltaTValue();

        return tmp<scalarField>(rDeltaT_());
    }
}

// ************************************************************************* //
