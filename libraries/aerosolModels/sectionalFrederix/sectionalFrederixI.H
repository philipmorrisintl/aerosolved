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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::List<Foam::scalar>
Foam::aerosolModels::sectionalFrederix::weights
(
    scalar& x1,
    scalar& x2,
    scalar& xs
)
{
    List<scalar> w(2, 0.0);

    // zeroth and third moment

    w[0] =  (xs-x2)/(x1-x2);
    w[1] = -(xs-x1)/(x1-x2);

    return w;
}

inline Foam::List<Foam::scalar>
Foam::aerosolModels::sectionalFrederix::weights
(
    scalar& x1,
    scalar& x2,
    scalar& x3,
    scalar& x4,
    scalar& xs
)
{
    List<scalar> w(4, 0.0);

    // zeroth, first, second and third moment

    w[0] = -((x2 - xs)*(x3 - xs)*(x4 - xs))/((x1 - x2)*(x1 - x3)*(x1 - x4));
    w[1] =  ((x1 - xs)*(x3 - xs)*(x4 - xs))/((x1 - x2)*(x2 - x3)*(x2 - x4));
    w[2] = -((x1 - xs)*(x2 - xs)*(x4 - xs))/((x1 - x3)*(x2 - x3)*(x3 - x4));
    w[3] =  ((x1 - xs)*(x2 - xs)*(x3 - xs))/((x1 - x4)*(x2 - x4)*(x3 - x4));

    return w;
}

inline void Foam::aerosolModels::sectionalFrederix::twoMoment
(
    Foam::label i,
    Foam::label j,
    Foam::scalar d,
    Foam::scalar G,
    Foam::scalar phi
)
{
    List<scalar> w = weights(x_[i], x_[i+1], d);

    M_[i][j]   += phi*w[0]*G;
    M_[i+1][j] += phi*w[1]*G;
}

inline void Foam::aerosolModels::sectionalFrederix::fourMoment
(
    Foam::label i,
    Foam::label j,
    Foam::scalar d,
    Foam::scalar G,
    Foam::scalar phi
)
{
    if ((i == 0) || (i == P_-2))
    {
        // Boundary points: revert to two-moment

        List<scalar> w = weights(x_[i], x_[i+1], d);

        M_[0][j] += phi*w[0]*G;
        M_[1][j] += phi*w[1]*G;
    }
    else
    {
        List<scalar> w = weights(x_[i-1], x_[i], x_[i+1], x_[i+2], d);

        label c = 0;

        for (label l = i-1; l <= min(i+2, P_-1); l++)
        {
            M_[l][j] += phi*w[c]*G;
            c++;
        }
    }
}

namespace Foam // Necessary for Doxygen
{

inline List<label> aerosolModels::sectionalFrederix::xLowerIndex
(
    List<scalar>& ds
) const
{
    List<label> k(ds.size(), -1);

    forAll(ds, i)
    {
        k[i] = xLowerIndex(ds[i]);
    }

    return k;
}

inline label aerosolModels::sectionalFrederix::xLowerIndex
(
    scalar ds
) const
{
    // Binary search algorithm. Return -1 if the diameter is too small and P_ if
    // the diameter is too large

    if (ds < x_[0])
    {
        return -1;
    }
    else if(ds > x_[P_-1])
    {
        return P_;
    }
    else
    {
        int imin = 0;
        int imax = P_-1;
        int imid = (imax+imin)/2;

        while (imin != imax-1)
        {
            if (x_[imid] >= ds)
            {
                imax = imid;
            }
            else
            {
                imin = imid;
            }

            imid = (imax+imin)/2;
        }

        return imin;
    }
}

} // End namespace Foam

// ************************************************************************* //
