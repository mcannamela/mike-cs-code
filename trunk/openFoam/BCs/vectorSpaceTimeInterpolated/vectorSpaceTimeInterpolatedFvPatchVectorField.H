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

Class
    vectorSpaceTimeInterpolatedFvPatchVectorField

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    vectorSpaceTimeInterpolatedFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef vectorSpaceTimeInterpolatedFvPatchVectorField_H
#define vectorSpaceTimeInterpolatedFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "linear3DLUT.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class vectorSpaceTimeInterpolatedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class vectorSpaceTimeInterpolatedFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        int t0Component;
		int t1Component;
		string LUTFile;
		linear_three_D_LUT LUT;
		//- Flow direction
        vector n_;


public:

    //- Runtime type information
    TypeName("vectorSpaceTimeInterpolated");


    // Constructors

        //- Construct from patch and internal field
        vectorSpaceTimeInterpolatedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        vectorSpaceTimeInterpolatedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given vectorSpaceTimeInterpolatedFvPatchVectorField
        //  onto a new patch
        vectorSpaceTimeInterpolatedFvPatchVectorField
        (
            const vectorSpaceTimeInterpolatedFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new vectorSpaceTimeInterpolatedFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        vectorSpaceTimeInterpolatedFvPatchVectorField
        (
            const vectorSpaceTimeInterpolatedFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new vectorSpaceTimeInterpolatedFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
