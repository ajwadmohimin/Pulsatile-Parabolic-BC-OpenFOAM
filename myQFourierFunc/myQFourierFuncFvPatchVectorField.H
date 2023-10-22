/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2023 MD. AJWAD MOHIMIN, CHITTAGONG UNIVERSITY OF ENGINEERING & TECHNOLOGY, BD
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

Class
    Foam::myQFourierFuncFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a Pulsatile Parabolic Velocity condition,
    calculated as:

        \f[
            U{t,r} = 2 * (Q{t} / (pi*R^2)) * (1 - (r/R)^2)
        \f]

    where
    \vartable
        U{t,r}  | Velocity at r radius on t time [m/s]
        Q{t}    | Flow rate on t time [m^3/s]
        R       | Maximum radius from center of the face [m]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd | Default
        Q | List of coefficients of flow curve decomposed Fourier series in complex form | yes |
        omega       | The angular frequency of the first harmonic of |
                    |a Fourier series of a periodic flow rate = 2*pi*T, T denotes Time period   | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       myQFourierFunc;
        Q       ( (3 0) (0 2) (2 4) ); // list of complex values
        omega   7.854;
        value      uniform (0 0 0); // optional dummy
    }
    \endverbatim

SourceFiles
    myQFourierFuncFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myQFourierFuncFvPatchVectorField_H
#define myQFourierFuncFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "complexField.H"
#include "complexVectorField.H"
#include "complex.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class myQFourierFuncFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class myQFourierFuncFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Fourier coefficients of Flow rate time series
        List<complex> Q_;

        //- Angular Frequency of first harmonic
        scalar omega_;

        //-Radius of the artery
        scalar R_;

    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("myQFourierFunc");


    // Constructors

        //- Construct from patch and internal field
        myQFourierFuncFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myQFourierFuncFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        myQFourierFuncFvPatchVectorField
        (
            const myQFourierFuncFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        myQFourierFuncFvPatchVectorField
        (
            const myQFourierFuncFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myQFourierFuncFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myQFourierFuncFvPatchVectorField
        (
            const myQFourierFuncFvPatchVectorField&,
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
                new myQFourierFuncFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            // //- Map (and resize as needed) from self given a mapping object
            // virtual void autoMap
            // (
            //     const fvPatchFieldMapper&
            // );
            //
            // //- Reverse map the given fvPatchField onto this fvPatchField
            // virtual void rmap
            // (
            //     const fvPatchVectorField&,
            //     const labelList&
            // );

        // Return Flow rate at t second
        complex Q_t(List<complex> Q_, scalar t)
        {
            complex Q_t_(0,0);
            complex ci(0,1);

            forAll(Q_, n)
            {
                Q_t_ += Q_[n] * Foam::exp(ci*n*omega_*t);
            }

            return Q_t_;
        }

        scalarField u_t(scalarField r_, scalar R_, scalar t)
        {
            complexField U_t_(r_.size());
            scalar pi(M_PI);

            forAll(r_,nr)
            {
                U_t_[nr] = 2*(Q_t(Q_, t)/(pi*R_*R_))*(1-Foam::pow((r_[nr]/R_),2));
            }
            return Re(U_t_);
        }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
