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

#include "icoFunction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::icoFunction<Specie>::icoFunction
(
    const Specie& sp,
    const PtrList<const Function1<scalar>>& functions,
    const scalarList& frac
)
:
    Specie(sp),
    functions_(0),
    frac_(frac)
{
    forAll(functions, i)
    {
        functions_.append(functions[i].clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::icoFunction<Specie>::icoFunction
(
    const icoFunction<Specie>& ip
)
:
    Specie(ip),
    functions_(0),
    frac_(ip.frac_)
{
    forAll(ip.functions_, i)
    {
        functions_.append(ip.functions_[i].clone().ptr());
    }
}


template<class Specie>
inline Foam::icoFunction<Specie>::icoFunction
(
    const word& name,
    const icoFunction<Specie>& ip
)
:
    Specie(name, ip),
    functions_(0),
    frac_(ip.frac_)
{
    forAll(ip.functions_, i)
    {
        functions_.append(ip.functions_[i].clone().ptr());
    }
}


template<class Specie>
inline Foam::autoPtr<Foam::icoFunction<Specie>>
Foam::icoFunction<Specie>::clone() const
{
    return autoPtr<icoFunction<Specie>>
    (
        new icoFunction<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::icoFunction<Specie>>
Foam::icoFunction<Specie>::New(const dictionary& dict)
{
    return autoPtr<icoFunction<Specie>>
    (
        new icoFunction<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    if (mag(sum(frac_)) < VSMALL)
    {
        return functions_[0].value(T);
    }
    else
    {
        scalar rho(0.0);

        forAll(functions_, i)
        {
            rho += frac_[i]*functions_[i].value(T);
        }

        return rho/sum(frac_);
    }
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::H
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoFunction<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::icoFunction<Specie>::operator=
(
    const icoFunction<Specie>& ip
)
{
    Specie::operator=(ip);

    functions_.clear();
    frac_.clear();

    forAll(ip.functions_, i)
    {
        functions_.append(ip.functions_[i].clone().ptr());
    }

    frac_ = ip.frac_;
}


template<class Specie>
inline void Foam::icoFunction<Specie>::operator+=
(
    const icoFunction<Specie>& ip
)
{
    forAll(ip.functions_, i)
    {
        functions_.append(ip.functions_[i].clone().ptr());
    }

    frac_.append(ip.frac_);
}


template<class Specie>
inline void Foam::icoFunction<Specie>::operator*=(const scalar s)
{
    frac_ = frac_*s;
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
Foam::icoFunction<Specie> Foam::operator+
(
    const icoFunction<Specie>& ip1,
    const icoFunction<Specie>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2)
    );

    PtrList<const Function1<scalar>> functions(0);

    forAll(ip1.functions_, i)
    {
        functions.append(ip1.functions_[i].clone().ptr());
    }

    forAll(ip2.functions_, i)
    {
        functions.append(ip2.functions_[i].clone().ptr());
    }

    scalarList frac(ip1.frac_);
    frac.append(ip2.frac_);

    return icoFunction<Specie>(sp, functions, frac);
}


template<class Specie>
Foam::icoFunction<Specie> Foam::operator*
(
    const scalar s,
    const icoFunction<Specie>& ip
)
{
    const scalarList frac(ip.frac_*s);

    return icoFunction<Specie>
    (
        s*static_cast<const Specie&>(ip),
        ip.functions_,
        frac
    );
}


// ************************************************************************* //
