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

#include "generalGaussianFixedGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

generalGaussianFixedGradientFvPatchScalarField::generalGaussianFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    maxVal_(1),
    offset_(0),
    n_(2),
    sx_(1),
    sy_(1),
    sz_(1),
    mux_(0),
    muy_(0),
    muz_(0)
{}


generalGaussianFixedGradientFvPatchScalarField::generalGaussianFixedGradientFvPatchScalarField
(
    const generalGaussianFixedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    maxVal_(ptf.maxVal_),
    offset_(ptf.offset_),
    n_(ptf.n_),
    sx_(ptf.sx_),
    sy_(ptf.sy_),
    sz_(ptf.sz_),
    mux_(ptf.mux_),
    muy_(ptf.muy_),
    muz_(ptf.muz_)
{}


generalGaussianFixedGradientFvPatchScalarField::generalGaussianFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    maxVal_(readScalar(dict.lookup("maxVal"))),
    offset_(readScalar(dict.lookup("offset"))),
    n_(readScalar(dict.lookup("n"))),
    sx_(readScalar(dict.lookup("sx"))),
    sy_(readScalar(dict.lookup("sy"))),
    sz_(readScalar(dict.lookup("sz"))),
    mux_(readScalar(dict.lookup("mux"))),
    muy_(readScalar(dict.lookup("muy"))),
    muz_(readScalar(dict.lookup("muz")))
{
    if (mag(sx_) < SMALL || mag(sy_) < SMALL || mag(sz_) < SMALL)
    {
        FatalErrorIn("generalGaussianFixedGradientFvPatchScalarField(dict)")
            << "small standard dev detected, that's gonna be a div by zero"
            << abort(FatalError);
    }


    evaluate();
}


generalGaussianFixedGradientFvPatchScalarField::generalGaussianFixedGradientFvPatchScalarField
(
    const generalGaussianFixedGradientFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    maxVal_(fcvpvf.maxVal_),
    offset_(fcvpvf.offset_),
    n_(fcvpvf.n_),
    sx_(fcvpvf.sx_),
    sy_(fcvpvf.sy_),
    sz_(fcvpvf.sz_),
    mux_(fcvpvf.mux_),
    muy_(fcvpvf.muy_),
    muz_(fcvpvf.muz_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void generalGaussianFixedGradientFvPatchScalarField::updateCoeffs()
{
        
    const vectorField& c = patch().Cf();
    vector mu(mux(), muy(), muz());
    
    //scalarField::operator=(offset()+maxVal()*exp(-pow(mag(c.component(0)-mux())/sx(),n())-pow(mag(c.component(1)-muy())/sy(),n())-pow(mag(c.component(2)-muz())/sz(),n())));
    
    scalarField::operator=(
        this->patchInternalField()+
                   (offset()+maxVal()*exp(-pow(
					         pow(
								  pow(mag(c.component(0)-mux())/sx(),2)
								 +pow(mag(c.component(1)-muy())/sy(),2)
								 +pow(mag(c.component(2)-muz())/sz(),2)
                             ,.5)
                             ,n()))
                           )/this->patch().deltaCoeffs()
        
        );
}


// Write
void generalGaussianFixedGradientFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    
    os.writeKeyword("maxVal")
        << maxVal_ << token::END_STATEMENT << nl;
    os.writeKeyword("offset")
        << offset_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sx")
        << sx_ << token::END_STATEMENT << nl;
    os.writeKeyword("sy")
        << sy_ << token::END_STATEMENT << nl;
    os.writeKeyword("sz")
        << sz_ << token::END_STATEMENT << nl;
        
    os.writeKeyword("mux")
        << mux_ << token::END_STATEMENT << nl;
    os.writeKeyword("muy")
        << muy_ << token::END_STATEMENT << nl;
    os.writeKeyword("muz")
        << muz_ << token::END_STATEMENT << nl;
        
    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, generalGaussianFixedGradientFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
