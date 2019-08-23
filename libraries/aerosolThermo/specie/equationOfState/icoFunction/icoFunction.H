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

/**

\file icoFunction.H
\brief equationOfState class which implements a custom Function1 as EoS

This class reads a Function1 object from the thermophysicalProperties file,
which will then act as the equation of state. This class can be used to
implement any equation of state. Normally, however, it is recommended to use
OpenFOAM's standard EoS classes, as the call to the Function1 object can be
quite expensive.

*/

#ifndef icoFunction_H
#define icoFunction_H

#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class icoFunction;

template<class Specie>
icoFunction<Specie> operator+
(
    const icoFunction<Specie>&,
    const icoFunction<Specie>&
);

template<class Specie>
icoFunction<Specie> operator*
(
    const scalar,
    const icoFunction<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const icoFunction<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class icoFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class icoFunction
:
    public Specie
{
    // Private data

        //- List of function objects
        PtrList<const Function1<scalar>> functions_;

        //- List of fractions
        scalarList frac_;


public:

    // Constructors

        //- Construct from components
        inline icoFunction
        (
            const Specie& sp,
            const PtrList<const Function1<scalar>>& functions,
            const scalarList& frac
        );

        //- Construct from dictionary
        icoFunction(const dictionary& dict);

        //- Construct as copy
        inline icoFunction(const icoFunction&);

        //- Construct as named copy
        inline icoFunction(const word& name, const icoFunction&);

        //- Construct and return a clone
        inline autoPtr<icoFunction> clone() const;

        // Selector from dictionary
        inline static autoPtr<icoFunction> New(const dictionary& dict);


    // Destructors

        ~icoFunction();


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoFunction<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const icoFunction&);
        inline void operator+=(const icoFunction&);
        inline void operator*=(const scalar);


    // Friend operators

        friend icoFunction operator+ <Specie>
        (
            const icoFunction&,
            const icoFunction&
        );

        friend icoFunction operator* <Specie>
        (
            const scalar s,
            const icoFunction&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const icoFunction&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeicoFunction                                                        \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    icoFunction<Specie>,                                                       \
    "icoFunction",                                                             \
    0                                                                          \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "icoFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
