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

\file pisoControl.H
\brief Library to control the PISO algorithm flow

This library controls how the PISO algorithm behaves. It reads its parameters
from the PISO dictionary in the fvSolution file. Some important functions are:

- The loop() function controls the main PISO corrector loop. It is usually
  called from a while argument. The loop() function returns true if convergence
  criteria are met. These criteria are:

  - The number of PISO correctors has exceeded the max number of PISO correctors
  - The first initial residual of the pressure equation is 1/relTol larger than
    the current initial residual of the pressure equation
  - The current initial residual of the pressure equation has dropped below the
    tolerance level.

  In order for the loop() function to be able to do this, the solver must
  register the initial residual of the pressure equation inside the pisoControl
  object. This is done using:

  \code{.unparsed}
  if (piso.corrNonOrtho() == 1)
  {
      piso.updateResidual();
  }
  \endcode

- The theta() function provides the value for \f$\theta\f$ in the so-called
  \f$\theta\f$-scheme.

- The corr() function provides the current corrector level (corresponds to the
  (k-1)th explicit corrector step as defined in Sec. 3.3 in \cite thesis.

- The correctNonOrthogonal() function offers an iterator for non-orthogonal
  correctors. A local `inner iteration' is performed with this function to make
  sure that explicit contributions inside an equation are updated sufficiently.
  This function is controlled by the nNonOrthogonalCorrectors directive in the
  PISO dict of the fvSolution file.

*/

#ifndef pisoControl_H
#define pisoControl_H

#include "solutionControl.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pisoControl Declaration
\*---------------------------------------------------------------------------*/

class pisoControl
:
    public solutionControl
{
protected:

    // Protected data

        // Solution controls

            //- Maximum number of PISO correctors
            label nCorrPISO_;

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- PISO corrector convergence criterum
            scalar residualCorrPISO_;

            //- PISO corrector relative tolerance
            scalar relTol_;

            //- PISO initial residual
            scalar initPISOresidual_;

            //- Switch to update material coefficients
            Switch updateCoeffsPISO_;

            //- Write residuals to file
            Switch writeResiduals_;

            //- File to write residuals to
            OFstream residualsOF_;

            //- Unscaled pressure equation residual
            scalar residual_;

            //- Theta scheme parameter
            scalar theta_;

            //- Info level memory
            label InfoLevel_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Disallow default bitwise copy construct
        pisoControl(const pisoControl&);

        //- Disallow default bitwise assignment
        void operator=(const pisoControl&);

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Write residual to log file
        void writeResidual();


public:


    // Static Data Members

        //- Run-time type information
        TypeName("pisoControl");


    // Constructors

        //- Construct from mesh
        pisoControl(fvMesh& mesh);


    //- Destructor
    virtual ~pisoControl();


    // Member Functions

        // Access

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- PISO corrector convergence criterum
            inline scalar residualCorrPISO() const;

            //- PISO corrector relative tolerance
            inline scalar relTol() const;

            //- Switch to update material coefficients
            inline Switch updateCoeffsPISO() const;

            //- Theta scheme parameter
            inline scalar theta() const;

        // Solution control

            //- PISO loop
            virtual bool loop();

            //- Check convergence criteria
            virtual bool converged();

            //- Store the residual
            void updateResidual();

            //- Get the current residual
            scalar getResidual();

            //- Silence the Info stream
            void silent();

            //- Reset the Info stream
            void verbal();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pisoControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
