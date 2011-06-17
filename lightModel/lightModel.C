/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "lightModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(lightModel, 0);
        defineRunTimeSelectionTable(lightModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::lightModel::lightModel(const volScalarField& intensity)
:
    IOdictionary
    (
        IOobject
        (
            "lightProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    light_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    absorptionEmission_(NULL),
    scatter_(NULL)
{}


Foam::light::lightModel::lightModel
(
    const word& type,
    const volScalarField& intensity
)
:
    IOdictionary
    (
        IOobject
        (
            "lightProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    light_(lookup("light")),
    coeffs_(subDict(type + "Coeffs")),
    solverFreq_(readLabel(lookup("solverFreq"))),
    absorptionEmission_(absorptionEmissionModel::New(*this, mesh_)),
    scatter_(scatterModel::New(*this, mesh_))
{
    solverFreq_ = max(1, solverFreq_);
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::light::lightModel::~lightModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::light::lightModel::read()
{
    if (regIOobject::read())
    {
        lookup("light") >> light_;
        coeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::light::lightModel::correct()
{
    if (!light_)
    {
        return;
    }

    if (time_.timeIndex() % solverFreq_ == 0)
    {
        calculate();
    }
}




// ************************************************************************* //
