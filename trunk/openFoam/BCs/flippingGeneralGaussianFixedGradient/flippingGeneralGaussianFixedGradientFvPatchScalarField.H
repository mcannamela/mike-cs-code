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
    flippingGeneralGaussianFixedGradientFvPatchScalarField

Description
    sets the normal component of the gradient equal to a fixed value in the shape of 
    a gaussian bump

SourceFiles
    flippingGeneralGaussianFixedGradientFvPatchScalarField.C

Author
    

\*---------------------------------------------------------------------------*/

#ifndef flippingGeneralGaussianFixedGradientFvPatchScalarField_H
#define flippingGeneralGaussianFixedGradientFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class flippingGeneralGaussianFixedGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class flippingGeneralGaussianFixedGradientFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        
        scalar maxVal_;
        scalar offset_;
        scalar n_;
        //- standard devs in all 3 directions
        scalar sx_;
        scalar sy_;
        scalar sz_;

        //- means in all 3 directions
        scalar mux_;
        scalar muy_;
        scalar muz_;
        
        scalar period_;
        scalar oldTime_;
        


public:

    //- Runtime type information
    TypeName("flippingGeneralGaussianFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        flippingGeneralGaussianFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flippingGeneralGaussianFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given flippingGeneralGaussianFixedGradientFvPatchScalarField
        //  onto a new patch
        flippingGeneralGaussianFixedGradientFvPatchScalarField
        (
            const flippingGeneralGaussianFixedGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flippingGeneralGaussianFixedGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flippingGeneralGaussianFixedGradientFvPatchScalarField
        (
            const flippingGeneralGaussianFixedGradientFvPatchScalarField&,
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
                new flippingGeneralGaussianFixedGradientFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return standard devs
        scalar& maxVal()
        {
            return maxVal_;
        }
        scalar& offset()
        {
            return offset_;
        }
        scalar& n()
        {
            return n_;
        }
        
        scalar& sx()
        {
            return sx_;
        }
        scalar& sy()
        {
            return sy_;
        }
        scalar& sz()
        {
            return sz_;
        }

        //- Return means devs
        scalar& mux()
        {
            return mux_;
        }
        scalar& muy()
        {
            return muy_;
        }
        scalar& muz()
        {
            return muz_;
        }
        vector mu()
        {
			vector mu(mux(), muy(), muz());
			return mu;
		}
        
       

       
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
