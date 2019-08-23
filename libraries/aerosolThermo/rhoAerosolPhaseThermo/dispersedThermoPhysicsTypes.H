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

#ifndef dispersedThermoPhysicsTypes_H
#define dispersedThermoPhysicsTypes_H

#include "specie.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "icoFunction.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on sensibleEnthalpy

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > constDispHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > dispHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                rhoConst<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleDispHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoFunction<specie>
            >,
            sensibleEnthalpy
        >
    > constIncompressibleFuncDispHThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleEnthalpy
        >
    > constIncompressiblePoly8DispHThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleEnthalpy
        >
    > adiabaticDispHThermoPhysics;

    // thermo physics types based on sensibleInternalEnergy

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > constDispEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > dispEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleDispEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoFunction<specie>
            >,
            sensibleInternalEnergy
        >
    > constIncompressibleFuncDispEThermoPhysics;

    typedef
    constTransport
    <
        species::thermo
        <
            hConstThermo
            <
                icoPolynomial<specie, 8>
            >,
            sensibleInternalEnergy
        >
    > constIncompressiblePoly8DispEThermoPhysics;

    typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                adiabaticPerfectFluid<specie>
            >,
            sensibleInternalEnergy
        >
    > adiabaticDispEThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
