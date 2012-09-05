/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "RCWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RCWallFvPatchScalarField::RCWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    R_(1),
    C_(1),
    T0_(300),
    k_(1),
    rhoC_(1000),
    oldt(-1)
{}


RCWallFvPatchScalarField::RCWallFvPatchScalarField
(
    const RCWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    R_(ptf.R_),
    C_(ptf.C_),
    T0_(ptf.T0_),
    k_(ptf.k_),
    rhoC_(ptf.rhoC_),
    oldt(ptf.oldt)
    
{}


RCWallFvPatchScalarField::RCWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    R_(readScalar(dict.lookup("R"))),
    C_(readScalar(dict.lookup("C"))),
    T0_(readScalar(dict.lookup("T0"))),
    k_(readScalar(dict.lookup("k"))),
    rhoC_(readScalar(dict.lookup("rhoC")))
{
    fvPatchField<scalar>::operator=
         (
             Field<scalar>("value", dict, p.size())
         );
    oldt = this->db().time().value();
    evaluate();
}


RCWallFvPatchScalarField::RCWallFvPatchScalarField
(
    const RCWallFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    R_(fcvpvf.R_),
    C_(fcvpvf.C_),
    T0_(fcvpvf.T0_),
    k_(fcvpvf.k_),
    rhoC_(fcvpvf.rhoC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void RCWallFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    double qin, Rair,Req,  Tw0, Tw, Tc, dt, t;
    dt  = this->db().time().deltaT().value();
    t = this->db().time().value();
    
    if (t==oldt)
		return;
    
    //this->db().time().value()
    word kName;
    double alphaConst = .0003;
    
    
    if (db().foundObject<volScalarField>("kappat"))
		{kName = "kappat";}
	else if (db().foundObject<volScalarField>("nut"))
	    {
			kName = "nut";
	    }
	else 
	    {kName = "T";}
	    
	Info << "kName is " <<kName <<endl;
	
	//const fvPatchField<scalar>& alpha = patch().lookupPatchField<volScalarField, scalar>(kName);
    
    const labelUList& faceCells = this->patch().faceCells();
    
    forAll((*this), facei)
    {
		//temperature and resistance of the adjoining cell
		Tc = this->internalField()[faceCells[facei]];
		if (kName == "T")
			Rair = 1.0/(this->patch().deltaCoeffs()[facei]*(k_+rhoC_*alphaConst));
		else
			Rair = 1.0/(this->patch().deltaCoeffs()[facei]*(k_+rhoC_*patch().lookupPatchField<volScalarField, scalar>(kName)[facei]));
		Req = Rair+R_*.5;
		
		//if (patch().lookupPatchField<volScalarField, scalar>(kName)[facei] > 1e-5)
		//	Info << patch().lookupPatchField<volScalarField, scalar>(kName)[facei] <<endl;
		
		(*this)[facei] += (2.0*dt/C_)*(Tc/Req+T0_*(1.0/R_-.5/Req));
		(*this)[facei] /= 1.0+(2.0*dt/C_)*(1.0/R_+.5/Req);
		
	}
	oldt = this->db().time().value();
}


// Write
void RCWallFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    
    os.writeKeyword("R")
        << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("C")
        << C_ << token::END_STATEMENT << nl;
    os.writeKeyword("T0")
        << T0_ << token::END_STATEMENT << nl;
    os.writeKeyword("k")
        << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoC")
        << rhoC_ << token::END_STATEMENT << nl;
        

    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, RCWallFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
