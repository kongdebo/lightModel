/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

Class
    Foam::semiTransSurfaceMixedFvPatchScalarField

Description
    light temperature specified

SourceFiles
    semiTransSurfaceMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef semiTransSurfaceMixedFvPatchScalarField_H
#define semiTransSurfaceMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{
/*---------------------------------------------------------------------------*\
        Class semiTransSurfaceMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class semiTransSurfaceMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- incident light angle distribution
        scalar  n1_,n2_;          // index of refraction 
        vector  d0_ ;         // incident light direction
        scalar  I0_;              // incident light intensity

public:

    //- Runtime type information
    TypeName("semiTransSurface");


    // Constructors

        //- Construct from patch and internal field
        semiTransSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        semiTransSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  semiTransSurfaceMixedFvPatchScalarField onto a new patch
        semiTransSurfaceMixedFvPatchScalarField
        (
            const semiTransSurfaceMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        semiTransSurfaceMixedFvPatchScalarField
        (
            const semiTransSurfaceMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new semiTransSurfaceMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        semiTransSurfaceMixedFvPatchScalarField
        (
            const semiTransSurfaceMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new semiTransSurfaceMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
