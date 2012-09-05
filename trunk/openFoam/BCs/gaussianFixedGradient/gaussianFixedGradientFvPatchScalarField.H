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
    gaussianFixedGradientFvPatchScalarField

Description
    sets the normal component of the gradient equal to a fixed value in the shape of 
    a gaussian bump

SourceFiles
    gaussianFixedGradientFvPatchScalarField.C

Author
    

\*---------------------------------------------------------------------------*/

#ifndef gaussianFixedGradientFvPatchScalarField_H
#define gaussianFixedGradientFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class gaussianFixedGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class gaussianFixedGradientFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        
        scalar maxVal_;
        scalar offset_;
        //- standard devs in all 3 directions
        scalar sx_;
        scalar sy_;
        scalar sz_;

        //- means in all 3 directions
        scalar mux_;
        scalar muy_;
        scalar muz_;
        
        
        


public:

    //- Runtime type information
    TypeName("gaussianFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        gaussianFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gaussianFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gaussianFixedGradientFvPatchScalarField
        //  onto a new patch
        gaussianFixedGradientFvPatchScalarField
        (
            const gaussianFixedGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gaussianFixedGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gaussianFixedGradientFvPatchScalarField
        (
            const gaussianFixedGradientFvPatchScalarField&,
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
                new gaussianFixedGradientFvPatchScalarField(*this, iF)
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
        
        //set means
        void set_mux(scalar mu)
        {
			mux_ = mu;
		}
		void set_muy(scalar mu)
        {
			muy_ = mu;
		}
		void set_muz(scalar mu)
        {
			muz_ = mu;
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
