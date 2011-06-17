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

Namespace
    Foam::light

Description
    Namespace for light modelling

Class
    Foam::light::lightModel

Description
    Top level model for light modelling

SourceFiles
    lightModel.C

\*---------------------------------------------------------------------------*/

#ifndef lightModel_H
#define lightModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{

// Forward declaration of classes
class absorptionEmissionModel;
class scatterModel;

/*---------------------------------------------------------------------------*\
                       Class lightModel Declaration
\*---------------------------------------------------------------------------*/

class lightModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;


        //- Model specific dictionary input parameters
        Switch light_;

        //- light model dictionary
        dictionary coeffs_;

        //- light solver frequency - number flow solver iterations per
        //  light solver iteration
        label solverFreq_;


        // References to the light sub-models

            //- Absorption/emission model
            autoPtr<absorptionEmissionModel> absorptionEmission_;

            //- Scatter model
            autoPtr<scatterModel> scatter_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        lightModel(const lightModel&);

        //- Disallow default bitwise assignment
        void operator=(const lightModel&);


public:

    //- Runtime type information
    TypeName("lightModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             lightModel,
             dictionary,
             (
                 const volScalarField& intensity
             ),
             (intensity)
         );


    // Constructors

        //- Null constructor
        lightModel(const volScalarField& intensity);

        //- Construct from components
        lightModel(const word& type, const volScalarField& intensity);


    // Selectors

         //- Return a reference to the selected light model
         static autoPtr<lightModel> New(const volScalarField& intensity);


    //- Destructor
    virtual ~lightModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve light equation(s)
            virtual void calculate() = 0;

            //- Read lightProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
  //          virtual tmp<volScalarField> Rp() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace light
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
