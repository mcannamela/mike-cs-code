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

#include "vectorSpaceTimeInterpolatedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vectorSpaceTimeInterpolatedFvPatchVectorField::vectorSpaceTimeInterpolatedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    t0Component(1),
    t1Component(2),
    n_(1, 0, 0)
{}

vectorSpaceTimeInterpolatedFvPatchVectorField::vectorSpaceTimeInterpolatedFvPatchVectorField
(
    const vectorSpaceTimeInterpolatedFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    t0Component(ptf.t0Component),
    t1Component(ptf.t1Component),
    LUTFile(ptf.LUTFile),
    LUT(ptf.LUT),
    n_(ptf.n_)
{}


vectorSpaceTimeInterpolatedFvPatchVectorField::vectorSpaceTimeInterpolatedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    t0Component(readScalar(dict.lookup("t0Component"))),
    t1Component(readScalar(dict.lookup("t1Component"))),
    n_(dict.lookup("n"))
{
    LUTFile = dict.lookup("LUTFile");
    LUT = linear_three_D_LUT(LUTFile);
   // LUT.printMe();

    evaluate();
}


vectorSpaceTimeInterpolatedFvPatchVectorField::vectorSpaceTimeInterpolatedFvPatchVectorField
(
    const vectorSpaceTimeInterpolatedFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    t0Component(fcvpvf.t0Component),
    t1Component(fcvpvf.t1Component),
    LUTFile(fcvpvf.LUTFile),
    LUT(fcvpvf.LUT),
    n_(fcvpvf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vectorSpaceTimeInterpolatedFvPatchVectorField::updateCoeffs()
{
    const vectorField& c = patch().Cf();
    
    
    forAll((*this), facei)
    {
		(*this)[facei] = (n_/mag(n_))*LUT(c[facei].component(t0Component), c[facei].component(t1Component), this->db().time().value());
	}

    
}


// Write
void vectorSpaceTimeInterpolatedFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("t0Component")
        << t0Component << token::END_STATEMENT << nl;
    os.writeKeyword("t1Component")
        << t1Component << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("LUTFile")
        << LUTFile << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, vectorSpaceTimeInterpolatedFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
