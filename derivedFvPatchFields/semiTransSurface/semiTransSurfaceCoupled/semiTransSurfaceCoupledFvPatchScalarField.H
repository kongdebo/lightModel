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

Class
    semiTransSurfaceCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used by the
    conjugate heat transfer solver.
    Both sides use a mix of zeroGradient and neighbourvalue.

    Example usage:
        myInterfacePatchName
        {
            type                semiTransSurfaceCoupled;
            neighbourFieldName  T;

        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

SourceFiles
    semiTransSurfaceCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef semiTransSurfaceCoupledFvPatchScalarField_H
#define semiTransSurfaceCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class semiTransSurfaceCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class semiTransSurfaceCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;
        

public:

    //- Runtime type information
    TypeName("semiTransSurfaceCoupled");


    // Constructors

        //- Construct from patch and internal field
        semiTransSurfaceCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        semiTransSurfaceCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  semiTransSurfaceCoupledFvPatchScalarField onto a new patch
        semiTransSurfaceCoupledFvPatchScalarField
        (
            const semiTransSurfaceCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new semiTransSurfaceCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        semiTransSurfaceCoupledFvPatchScalarField
        (
            const semiTransSurfaceCoupledFvPatchScalarField&,
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
                new semiTransSurfaceCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
