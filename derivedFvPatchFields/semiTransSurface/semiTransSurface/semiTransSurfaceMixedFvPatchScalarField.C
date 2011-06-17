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

#include "semiTransSurfaceMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "lightDOM.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::light::semiTransSurfaceMixedFvPatchScalarField::
semiTransSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    n1_(0.0),
    n2_(0.0),
    d0_(vector::zero),
    I0_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::light::semiTransSurfaceMixedFvPatchScalarField::
semiTransSurfaceMixedFvPatchScalarField
(
    const semiTransSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    n1_(ptf.n1_),
    n2_(ptf.n2_),
    d0_(ptf.d0_),
    I0_(ptf.I0_)
{}


Foam::light::semiTransSurfaceMixedFvPatchScalarField::
semiTransSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    n1_(readScalar(dict.lookup("n1"))), //scalarField("n1", dict),
	n2_(readScalar(dict.lookup("n2"))), // scalarField("n2", dict);
  	d0_(dict.lookup("d0")), // scalarField("d0", dict);
    I0_(readScalar(dict.lookup("I0")))   //scalarField("I0", dict) ;
{


    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        
        
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
              
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
        
    }
    
//    Info << "\n n1  \t" << n1_ << endl;

    
}


Foam::light::semiTransSurfaceMixedFvPatchScalarField::
semiTransSurfaceMixedFvPatchScalarField
(
    const semiTransSurfaceMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    n1_(ptf.n1_),
    n2_(ptf.n2_),
    d0_(ptf.d0_),
    I0_(ptf.I0_)
{}


Foam::light::semiTransSurfaceMixedFvPatchScalarField::
semiTransSurfaceMixedFvPatchScalarField
(
    const semiTransSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    n1_(ptf.n1_),
    n2_(ptf.n2_),
    d0_(ptf.d0_),
    I0_(ptf.I0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::light::semiTransSurfaceMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

      
    scalarField& Iw = *this;
    
    const lightModel& light = db().lookupObject<lightModel>("lightProperties");

    const lightDOM& dom(refCast<const lightDOM>(light));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);
    
    if (dom.nLambda() != 1)
    {	FatalErrorIn("semiTransSurfaceMixedFvPatchScalarField::updateCoeffs") 
          << "absorption model" << nl << exit(FatalError);  }
    
    vectorField n = patch().Sf()/patch().magSf();
    
//   const label patchI = patch().index();
//   lightIntensityRay& ray = const_cast<lightIntensityRay&>(dom.IRay(rayId));
//    const scalarField& IFace =  dom.IRay(rayId).ILambda(lambdaId).boundaryField()[patchI];  
//   const vector& dAve = dom.IRay(rayId).dAve(); 
//    ray.Qr().boundaryField()[patchI] += Iw*(n & ray.dAve());    

    const vector& d = dom.IRay(rayId).dAve();
   

    scalar sinTheta1 ,sinTheta2,sinTheta0;
    scalar  r1,r2,R;
    scalar cosTheta0,cosTheta1,cosTheta2;

/*
   
    Info << "\n n1  \t" << n1_ << endl;
	Info << "\n n2  \t" << n2_ << endl;
	
  	Info << "\n d0  \t" << d0_ << endl;
    Info << "\n I0  \t" << I0_ << endl;
    */
    
    forAll(Iw, faceI)
    {

        cosTheta0 = -n[faceI] & d0_ ;  
        cosTheta2 = -n[faceI] & d  ; // / mag(d);                //    /mag(n[faceI])
        sinTheta2 = mag(-n[faceI] ^ d )  ;  //sinTheta2 
       sinTheta0  = mag(-n[faceI] ^ d0_); //  sinTheta0 

         
  
        if ( cosTheta2 > 0.0  &&  cosTheta0 > 0.0)   // direction out of the wall
        {
            sinTheta2 = mag(-n[faceI] ^ d)   ; //  / mag(d);       //  /mag(n[faceI])
            sinTheta1 = n2_ * sinTheta2 / n1_ ; 
            
            if(sinTheta1 <= 1 )
            {

            cosTheta1 = Foam::cos(Foam::asin(sinTheta1));
            
      //      if ( mag(Foam::asin(sinTheta0) - Foam::asin(sinTheta1))*180/Foam::mathematicalConstant::pi < 5 )
            {
				r1 = (n1_*cosTheta1 - n2_*cosTheta2 )
				  /  (n1_*cosTheta1 + n2_*cosTheta2 );
				r2 = (n2_*cosTheta1 - n1_*cosTheta2 )
				  /  (n1_*cosTheta2 + n2_*cosTheta1 );
            
				R = 0.5*(r1*r1 + r2*r2);
                        
				refValue()[faceI] = I0_ * Foam::cos(Foam::asin(sinTheta0) - Foam::asin(sinTheta1)) * (1 - R) ;         
			} 
	//		else
			{
//				refValue()[faceI] = 0.0 ;        
			}
            
            refGrad()[faceI] = 0.0;  //not used
            valueFraction()[faceI] = 1.0;
		    }

        }
        else                         // direction into the wall   
        {
            
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used
                   
        }
        

    }
    
    
    mixedFvPatchScalarField::updateCoeffs(); 
    


}


void Foam::light::semiTransSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("n1") << n1_ << token::END_STATEMENT << nl;
    os.writeKeyword("n2") << n2_ << token::END_STATEMENT << nl;
    os.writeKeyword("d0") << d0_ << token::END_STATEMENT << nl;
    os.writeKeyword("I0") << I0_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace light
{
    makePatchTypeField
    (
        fvPatchScalarField,
        semiTransSurfaceMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
