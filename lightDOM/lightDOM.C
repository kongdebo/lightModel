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

\*---------------------------------------------------------------------------*/

#include "lightDOM.H"
#include "addToRunTimeSelectionTable.H"


#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace light
    {
        defineTypeNameAndDebug(lightDOM, 0);

        addToRunTimeSelectionTable
        (
            lightModel,
            lightDOM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::lightDOM::lightDOM(const volScalarField& intensity)
:
    lightModel(typeName, intensity),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
     //   dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
         dimensionedScalar("G", dimLuminousIntensity, 0.0)
    ),
    
  /*  Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ), */  
    
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
 
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50))
{
    if (mesh_.nSolutionD() == 3)    //3D
    {
        nRay_ = 4*nPhi_*nTheta_;
        IRay_.setSize(nRay_);

        scalar deltaPhi   =  pi /(2.0*nPhi_);
        scalar deltaTheta =  pi  /nTheta_;

        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new lightIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_
                    )
                );
                i++;
            }
        }
    }
    else
    {
        if (mesh_.nSolutionD() == 2)    //2D (X & Y)
        {

            scalar thetai =      piByTwo;
            scalar deltaTheta =   pi;
            nRay_ = 4*nPhi_;
            IRay_.setSize(nRay_);
            scalar deltaPhi =    pi /(2.0*nPhi_);

            label i = 0;
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new lightIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_
                    )
                );
                i++;
            }
        }
        else    //1D (X)
        {

            scalar thetai =       piByTwo;
            scalar deltaTheta =   pi;
            nRay_ = 2;
            IRay_.setSize(nRay_);
            scalar deltaPhi =  pi;

            label i = 0;
            for (label m = 1; m <= 2; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new lightIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_
                    )
                );
                i++;
            }

        }
    }


    // Construct absorptionEmission field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "lightDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;
    forAll (IRay_, i)
    {
        Info<< '\t' << IRay_[i].I().name()
            << '\t' << IRay_[i].dAve() << nl;
    }
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::light::lightDOM::~lightDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::light::lightDOM::read()
{
    if (lightModel::read())
    {
//      Only reading solution parameters - not changing ray geometry

        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::light::lightDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    scalar maxResidual = 0.0;
    label radIter = 0;
    do
    {
        radIter++;
        forAll(IRay_, rayI)
        {
            maxResidual = 0.0;
            scalar maxBandResidual = IRay_[rayI].correct();
            maxResidual = max(maxBandResidual, maxResidual);
        }

        Info << "light solver iter: " << radIter << endl;

    } while(maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}



void Foam::light::lightDOM::updateG()
{
 //   G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
     G_ = dimensionedScalar("zero",dimLuminousIntensity, 0.0);
     
  //  Qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        
        //Qr_ += IRay_[rayI].Qr();
//        Qr_.boundaryField() += IRay_[rayI].Qr().boundaryField();
    }
}


void Foam::light::lightDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // assuming name is in the form: CHARS_rayId_lambdaId
    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    rayId = readLabel(IStringStream(name.substr(i1+1, i2-1))());
    lambdaId = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());
}


// ************************************************************************* //
