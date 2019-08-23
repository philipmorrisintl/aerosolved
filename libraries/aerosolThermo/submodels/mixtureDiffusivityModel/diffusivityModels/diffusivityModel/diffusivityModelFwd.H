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

#ifndef diffusivityModelFws_H
#define diffusivityModelFws_H

#include "diffusivityModel.H"
#include "vector.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef diffusivityModel<label> labeldiffusivityModel;
    typedef diffusivityModel<scalar> scalardiffusivityModel;
    typedef diffusivityModel<vector> vectordiffusivityModel;
    typedef diffusivityModel<symmTensor> symmTensordiffusivityModel;
    typedef diffusivityModel<sphericalTensor> sphericalTensordiffusivityModel;
    typedef diffusivityModel<tensor> tensordiffusivityModel;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
