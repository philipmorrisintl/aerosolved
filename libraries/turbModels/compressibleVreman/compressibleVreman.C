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

#include "compressibleVreman.H"
#include "makeTurbulenceTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbModels
{
    makeTurbModelTypes(compressibleVreman, turbModel);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbModels::compressibleVreman::compressibleVreman
(
    const fvMesh& mesh,
    fluidThermo& thermo
)
:
    turbModel(mesh, thermo)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbModels::compressibleVreman::~compressibleVreman()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbModels::compressibleVreman::update()
{
    volTensorField gradU = fvc::grad(U());

    volTensorField Beta  = (gradU.T() & gradU);
    volScalarField abeta = tr(Beta);
    volScalarField bbeta = Beta.component(tensor::XX)*Beta.component(tensor::YY)
                         - Beta.component(tensor::XY)*Beta.component(tensor::XY)
                         + Beta.component(tensor::XX)*Beta.component(tensor::ZZ)
                         - Beta.component(tensor::XZ)*Beta.component(tensor::XZ)
                         + Beta.component(tensor::YY)*Beta.component(tensor::ZZ)
                         - Beta.component(tensor::YZ)*Beta.component(tensor::YZ);

    if (Foam::min(bbeta) < dimensionedScalar("zero",bbeta.dimensions(),0.0))
    {
        Info<< "Warning: bbeta = " << Foam::min(bbeta).value() << endl;
    }

    volScalarField bbetasqrt = Foam::sqrt( Foam::sqrt(bbeta*bbeta));
    muTurb_ = (5.0/2.0)*rho()*cs_*cs_*delta()*delta()*(bbetasqrt)/Foam::sqrt(abeta + dimensionedScalar("VSMALLDS",abeta.dimensions(),VSMALL));
    muTurb_.correctBoundaryConditions();

    kTurb_ = thermo_.CpEff()*muTurb()/Prt();
    kTurb_.correctBoundaryConditions();

}


bool Foam::turbModels::compressibleVreman::read()
{
    delta_.internalField() = deltaCoeff()*pow(mesh_.V(), 1.0/3.0);
    delta_.correctBoundaryConditions();

    if (regIOobject::read())
    {
       coeffs_.lookup("cs") >> cs_;
       coeffs_.lookup("c1") >> c1_;
       coeffs_.lookup("c2") >> c2_;
       coeffs_.lookup("Prt") >> Prt_;
       coeffs_.lookup("deltaCoeff") >> deltaCoeff_;

       Info << "LES: Cs= " << cs_ << " c1= " << c1_ << " c2= "
            << c2_ <<" Prt= " << Prt_ << " delta= " << deltaCoeff_ << endl;

       return true;
    }
    else
    {
        return false;
    }

    return true;
}

Foam::tmp<Foam::volScalarField> Foam::turbModels::compressibleVreman::K(volVectorField& U )
{
    volTensorField gradU = fvc::grad(U);

    return
    (
        2.0*c2()*muTurb()*tr(symm(gradU))/3.0 - c1()*sqrt(2.0)*muTurb()*mag(symm(gradU))/3.0
    );
}


// ************************************************************************* //
