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

\file coeffs.H
\brief Update coefficients \f$C_p\f$, \f$\kappa\f$ and \f$\mu\f$

Update heat capacity at constant pressure \f$C_p\f$, thermal conductivity
\f$\kappa\f$, viscosity \f$\mu\f$ and the turbulence terms.

*/

thermo.updateCpEff();
conductivity.update();
viscosity.update();
turbulence.update();
