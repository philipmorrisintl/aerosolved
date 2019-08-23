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

inline label sectionalDistribution::findNearest
(
    const scalar& s,
    const bool soft
) const
{
    if (s<yMin())
    {
        return soft ? 0 : -1;
    }
    else if (s>=yMax())
    {
        return soft ? size()-1 : -1;
    }
    else
    {
        return min(max(search(s),0),size()-1);
    }
}

inline label sectionalDistribution::findLower
(
    const scalar& s,
    const bool soft
) const
{
    if (s<xMin())
    {
        return soft ? 0 : -1;
    }
    else if (s>=xMax())
    {
        return soft ? size()-2 : -1;
    }
    else
    {
        return min(max(searchLower(s),0),size()-2);
    }
}

inline label sectionalDistribution::findUpper
(
    const scalar& s,
    const bool soft
) const
{
    if (s<xMin())
    {
        return soft ? 1 : -1;
    }
    else if (s>=xMax())
    {
        return soft ? size()-1 : -1;
    }
    else
    {
        return min(max(searchLower(s)+1,1),size()-1);
    }
}

inline scalar sectionalDistribution::median
(
    const label& celli,
    const scalar& p
) const
{
    scalarList pM(size(), 0.0);

    pM[0] =
        pow(sections()[0].x(), p/3.0)
      * sections()[0].M().field()[celli];

    for (label i = 1; i < size(); i++)
    {
        pM[i] =
            pM[i-1]
          + pow(sections()[i].x(), p/3.0)
          * sections()[i].M().field()[celli];
    }

    pM = pM/max(pM[size()-1], VSMALL);

    if (pM[size()-1] == 1.0)
    {
        label j = Foam::findLower(pM, 0.5);

        if (j == -1)
        {
            return 0.0;
        }
        else if (j == size())
        {
            return sections()[size()-1].x();
        }
        else
        {
            return
                sections()[j].x()
              + (0.5-pM[j])/(pM[j+1]-pM[j])
              * (sections()[j+1].x()-sections()[j].x());
        }
    }
    else
    {
        return 0.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
