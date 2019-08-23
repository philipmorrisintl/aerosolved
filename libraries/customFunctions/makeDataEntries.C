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

#include "exponential.H"

#include "NSRDS0.H"
#include "NSRDS1.H"
#include "NSRDS2.H"
#include "NSRDS3.H"
#include "NSRDS4.H"
#include "NSRDS5.H"
#include "NSRDS6.H"
#include "NSRDS7.H"

#include "VDI1.H"
#include "VDI2.H"
#include "VDI3.H"
#include "VDI4.H"
#include "VDI5.H"
#include "VDI6.H"
#include "VDI7.H"
#include "VDI8.H"
#include "VDI9.H"
#include "VDI10.H"

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1s(Type)                                                   \
    makeFunction1Type(exponential, Type);

namespace Foam
{
    makeFunction1s(scalar);
    makeFunction1s(vector);
    makeFunction1s(sphericalTensor);
    makeFunction1s(symmTensor);
    makeFunction1s(tensor);

    makeFunction1Type(NSRDS0, scalar);
    makeFunction1Type(NSRDS1, scalar);
    makeFunction1Type(NSRDS2, scalar);
    makeFunction1Type(NSRDS3, scalar);
    makeFunction1Type(NSRDS4, scalar);
    makeFunction1Type(NSRDS5, scalar);
    makeFunction1Type(NSRDS6, scalar);
    makeFunction1Type(NSRDS7, scalar);

    makeFunction1Type(VDI1, scalar);
    makeFunction1Type(VDI2, scalar);
    makeFunction1Type(VDI3, scalar);
    makeFunction1Type(VDI4, scalar);
    makeFunction1Type(VDI5, scalar);
    makeFunction1Type(VDI6, scalar);
    makeFunction1Type(VDI7, scalar);
    makeFunction1Type(VDI8, scalar);
    makeFunction1Type(VDI9, scalar);
    makeFunction1Type(VDI10, scalar);
}

// ************************************************************************* //
