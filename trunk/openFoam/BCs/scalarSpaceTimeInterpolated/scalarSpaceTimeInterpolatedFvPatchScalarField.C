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

#include "scalarSpaceTimeInterpolatedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scalarSpaceTimeInterpolatedFvPatchScalarField::scalarSpaceTimeInterpolatedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    t0Component(1),
    t1Component(2)
{}


scalarSpaceTimeInterpolatedFvPatchScalarField::scalarSpaceTimeInterpolatedFvPatchScalarField
(
    const scalarSpaceTimeInterpolatedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    t0Component(ptf.t0Component),
    t1Component(ptf.t1Component),
    LUTFile(ptf.LUTFile),
    LUT(ptf.LUT)
{}


scalarSpaceTimeInterpolatedFvPatchScalarField::scalarSpaceTimeInterpolatedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    t0Component(readScalar(dict.lookup("t0Component"))),
    t1Component(readScalar(dict.lookup("t1Component")))
    
{
    LUTFile = dict.lookup("LUTFile");
    LUT = linear_three_D_LUT(LUTFile);
    //LUT.printMe();

    evaluate();
}


scalarSpaceTimeInterpolatedFvPatchScalarField::scalarSpaceTimeInterpolatedFvPatchScalarField
(
    const scalarSpaceTimeInterpolatedFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    t0Component(fcvpvf.t0Component),
    t1Component(fcvpvf.t1Component),
    LUTFile(fcvpvf.LUTFile),
    LUT(fcvpvf.LUT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void scalarSpaceTimeInterpolatedFvPatchScalarField::updateCoeffs()
{
    
    const vectorField& c = patch().Cf();
    forAll((*this), facei)
    {
		(*this)[facei] = LUT(c[facei].component(t0Component), c[facei].component(t1Component), this->db().time().value());
		
	}

}


// Write
void scalarSpaceTimeInterpolatedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("t0Component")
        << t0Component << token::END_STATEMENT << nl;
    os.writeKeyword("t1Component")
        << t1Component << token::END_STATEMENT << nl;
    os.writeKeyword("LUTFile")
        << LUTFile << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, scalarSpaceTimeInterpolatedFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
