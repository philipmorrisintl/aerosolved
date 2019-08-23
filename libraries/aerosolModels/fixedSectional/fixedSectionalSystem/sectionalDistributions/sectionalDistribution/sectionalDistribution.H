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

\file sectionalDistribution.H
\brief Base class for the distribution of sections

This class contains a list of sections, and determines the representative sizes
of those sections using a specified methods.

*/

#ifndef sectionalDistribution_H
#define sectionalDistribution_H

#include "aerosolSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "section.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sectionalDistribution Declaration
\*---------------------------------------------------------------------------*/

class sectionalDistribution
:
    public PtrList<section>,
    public aerosolSubModelBase
{
private:

    // Private Data

        //- Dimension of the size scalar
        dimensionSet sizeDims_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sectionalDistribution(const sectionalDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const sectionalDistribution&);


protected:

    // Protected Data

        //- Representative size list
        scalarList x_;

        //- Edge list
        scalarList y_;

public:

    //- Runtime type information
    TypeName("sectionalDistribution");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             sectionalDistribution,
             dictionary,
             (
                aerosolModel& aerosol,
                const dictionary& dict,
                const dimensionSet& dimensions
             ),
             (aerosol, dict, dimensions)
         );

    // Constructors

        //- Construct from type name, dictionary and aerosol model
        sectionalDistribution
        (
            const word& modelType,
            aerosolModel& aerosol,
            const dictionary& dict,
            const dimensionSet& dimensions
        );


    // Selectors

        //- Return a reference to the selected sectional distribution
        static autoPtr<sectionalDistribution> New
        (
            aerosolModel& aerosol,
            const dictionary& dict,
            const dimensionSet& dimensions
        );


    //- Destructor
    virtual ~sectionalDistribution();

    //- Access

        //- Const access to the size scalar dimensions
        inline const dimensionSet& sizeDimensions() const
        {
            return sizeDims_;
        }

        //- Access to the list of sections
        inline PtrList<section>& sections()
        {
            return *this;
        }

        //- Const access to the list of sections
        inline const PtrList<section>& sections() const
        {
            return *this;
        }

    // Member Functions

        //- Name of a given section
        word sectionName(const label i) const;

        //- List of representative sizes
        inline const scalarList& x() const
        {
            return x_;
        }

        //- List of section edges
        inline const scalarList& y() const
        {
            return y_;
        }

        //- Minimum of distribution size range
        inline scalar yMin() const
        {
            return y_[0];
        }

        //- Maximum of distribution size range
        inline scalar yMax() const
        {
            return y_[y_.size()-1];
        }

        //- Minimum of representative sizes
        inline scalar xMin() const
        {
            return x_[0];
        }

        //- Maximum of representative sizes
        inline scalar xMax() const
        {
            return x_[x_.size()-1];
        }

        //- Search the section number in which given size is located
        virtual label search(const scalar& s) const = 0;

        //- Search the section number of which the representative size is
        //  immediately left
        virtual label searchLower(const scalar& s) const = 0;

        //- Find the section number in which the given size is located
        inline label findNearest
        (
            const scalar& s,
            const bool soft = false
        ) const;

        //- Find the section number of which the representative size is
        //  immediately left
        inline label findLower
        (
            const scalar& s,
            const bool soft = false
        ) const;

        //- Find the section number of which the representative size is
        //  immediately right
        inline label findUpper
        (
            const scalar& s,
            const bool soft = false
        ) const;

        //- Compute the p moment median size
        inline scalar median
        (
            const label& cellI,
            const scalar& p
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sectionalDistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
