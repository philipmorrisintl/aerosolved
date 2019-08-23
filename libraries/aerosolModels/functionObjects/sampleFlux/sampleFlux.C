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

#include "sampleFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sampleFlux, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sampleFlux::sampleFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name, dict),
    patchNames_
    (
        dict.lookupOrDefault<wordList>("patches", wordList(0))
    ),
    faceZoneSetNames_
    (
        dict.lookupOrDefault<wordList>("faceZoneSets", wordList(0))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sampleFlux::~sampleFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sampleFlux::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && logFiles::read(dict))
    {
        wordList all(0);

        forAll(patchNames_, patchi)
        {
            forAll(mesh_.boundaryMesh(), patchj)
            {
                if (mesh_.boundaryMesh()[patchj].name() == patchNames_[patchi])
                {
                    all.append
                    (
                        IOobject::groupName("patch", patchNames_[patchi])
                    );

                    patchIDs_.append(patchj);

                    break;
                }
            }
        }

        forAll(faceZoneSetNames_, seti)
        {
            forAll(mesh_.faceZones(), setj)
            {
                if (mesh_.faceZones()[setj].name() == faceZoneSetNames_[seti])
                {
                    all.append
                    (
                        IOobject::groupName
                        (
                            "faceZoneSet",
                            faceZoneSetNames_[seti]
                        )
                    );

                    faceZoneSetIDs_.append(setj);

                    break;
                }
            }
        }

        logFiles::resetNames(all);

        Info<< type() << " " << name() << ": ";

        if (writeToFile() && names().size())
        {
            Info<< "applying to patches and faceZoneSets:" << nl;

            forAll(names(), i)
            {
                Info<< "    " << names()[i] << nl;
                writeFileHeader(files(i));
            }

            Info<< endl;
        }
        else
        {
            Info<< "no patches or faceZoneSets to be processed" << nl << endl;
        }

        return true;
    }

    return true;
}


bool Foam::functionObjects::sampleFlux::execute()
{
    return true;
}


bool Foam::functionObjects::sampleFlux::write()
{
    if (Pstream::master())
    {
        forAll(names(), filei)
        {
            writeTime(files(filei));
        }
    }

    forAll(fluxFieldNames_, fluxi)
    {
        const word fluxFieldName(fluxFieldNames_[fluxi]);

        if (mesh_.foundObject<surfaceScalarField>(fluxFieldName))
        {
            const surfaceScalarField& phi =
                mesh_.lookupObject<surfaceScalarField>(fluxFieldName);

            forAll(patchIDs_, patchi)
            {
                const label patchj(patchIDs_[patchi]);
                const label filei(patchi);

                const scalarField& phip = phi.boundaryField()[patchj];

                const scalar flux(gSum(phip));

                if (Pstream::master())
                {
                    files(filei) << token::TAB << flux;
                }
            }

            forAll(faceZoneSetIDs_, seti)
            {
                const label setj(faceZoneSetIDs_[seti]);
                const label filei(patchIDs_.size()+seti);

                scalar flux(0.0);

                forAll(mesh_.faceZones()[setj], facei)
                {
                    const label facej(mesh_.faceZones()[setj][facei]);

                    flux += phi[facej];
                }

                reduce(flux, sumOp<scalar>());

                if (Pstream::master())
                {
                    files(filei) << token::TAB << flux;
                }
            }
        }
    }

    if (Pstream::master())
    {
        forAll(names(), filei)
        {
            files(filei) << endl;
        }
    }

    return true;
}


// ************************************************************************* //
