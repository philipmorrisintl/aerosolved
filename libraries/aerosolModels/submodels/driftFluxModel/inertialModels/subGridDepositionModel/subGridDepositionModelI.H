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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool subGridDepositionModel::collision() const
{
    if (mag(u_) < tolerance_)
    {
        if
        (
            g_ < 0.0
         || (mag(g_) < tolerance_ && 1.0+v0_ < r_)
         || (
                mag(g_) >= tolerance_
             && (g_-v0_)/g_ > 1.0
             && x(log((g_-v0_)/g_)) < r_
            )
        )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        const scalar D(sqrt(max(4.0*u_+1.0,0.0)));

        const scalar l1(-D-1);
        const scalar l2(D-1);

        const scalar k(g_/u_);
        const scalar q(1.0+k);

        if (mag(D) < tolerance_)
        {
            if
            (
                -k < r_
             || (
                    mag(v0_+q/2.0) > tolerance_
                 && 2.0*v0_/(v0_+q/2.0) > 0.0
                 && x(2.0*v0_/(v0_+q/2.0)) < r_
                )
            )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            const scalar c
            (
                l1*(l2*q-2.0*v0_)
              / stabilise(l2*(l1*q-2.0*v0_), tolerance_)
            );

            if
            (
                (l2 < 0.0 && -k < r_)
             || (l2 > 0.0 && (l1*q-2.0*v0_) > 0.0)
             || (
                    mag(l1*q-2.0*v0_) > tolerance_
                 && c > 1.0
                 && x(2.0*log(1.0/c)/(l1-l2)) < r_
                )
            )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}

inline scalar subGridDepositionModel::x(const scalar t) const
{
    if (mag(u_) < tolerance_)
    {
        return (g_-v0_)*(exp(-t)-1.0)+1.0+g_*t;
    }
    else
    {
        const scalar D(sqrt(max(4.0*u_+1.0,0.0)));

        const scalar l1(-D-1);
        const scalar l2(D-1);

        const scalar k(g_/u_);
        const scalar q(1.0+k);

        if (mag(D) < tolerance_)
        {
            return (q+(v0_+q/2.0)*t)*exp(-t/2.0)-k;
        }
        else
        {
            return
              - (l2*q-2.0*v0_)/(l1-l2)*exp(min(t*l1/2.0,1e2))
              + (l1*q-2.0*v0_)/(l1-l2)*exp(min(t*l2/2.0,1e2))
              - k;
        }
    }
}

inline scalar subGridDepositionModel::v(const scalar t) const
{
    if (mag(u_) < tolerance_)
    {
        return -(g_-v0_)*exp(-t)+g_;
    }
    else
    {
        const scalar D(sqrt(max(4.0*u_+1.0,0.0)));

        const scalar l1(-D-1);
        const scalar l2(D-1);

        const scalar k(g_/u_);
        const scalar q(1.0+k);

        if (mag(D) < tolerance_)
        {
            return -0.5*(q+(v0_+q/2.0)*(t-2.0))*exp(-t/2.0);
        }
        else
        {
            return
              - 0.5*l1*(l2*q-2.0*v0_)/(l1-l2)*exp(min(t*l1/2.0,1e2))
              + 0.5*l2*(l1*q-2.0*v0_)/(l1-l2)*exp(min(t*l2/2.0,1e2));
        }
    }
}

inline scalar subGridDepositionModel::t() const
{
    const scalar t1(1.0/max(tolerance_, mag(u_)));
    const scalar t2(1.0/max(tolerance_, mag(v0_)));
    const scalar t3(1.0/max(tolerance_, mag(g_)));

    scalar t0(min(t1, min(t2, t3)));
    scalar t(0.0);

    for (label i = 0; i <= maxIter_; i++)
    {
        t = max(min(t0 - (x(t0)-r_)/v(t0), 1e2), -1e2);

        if (mag((t-t0)/t) < tolerance_)
        {
            break;
        }
        else if (i == maxIter_)
        {
            WarningInFunction
                << "Soluton did not converge"
                << endl;
        }
        else
        {
            t0 = t;
        }
    }

    return t;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
