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

#include "error.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::light::absorptionEmissionModel>
Foam::light::absorptionEmissionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word absorptionEmissionModelType(dict.lookup("absorptionEmissionModel"));

    Info<< "Selecting absorptionEmissionModel "
        << absorptionEmissionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(absorptionEmissionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "absorptionEmissionModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown absorptionEmissionModelType type "
            << absorptionEmissionModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid absorptionEmissionModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<absorptionEmissionModel>(cstrIter()(dict, mesh));
}


// ************************************************************************* //
