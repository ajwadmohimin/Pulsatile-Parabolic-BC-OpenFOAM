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

\*---------------------------------------------------------------------------*/

#include "myQFourierFuncFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::myQFourierFuncFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myQFourierFuncFvPatchVectorField::
myQFourierFuncFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Q_(2,complex {0,0}),
    omega_(0.0),
    R_(0.0)
{}


Foam::myQFourierFuncFvPatchVectorField::
myQFourierFuncFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Q_(dict.lookup("Q")),
    omega_(dict.get<scalar>("omega"))
{
    Info << "Using the myQFourierFunc boundary condition" << endl;

}


Foam::myQFourierFuncFvPatchVectorField::
myQFourierFuncFvPatchVectorField
(
    const myQFourierFuncFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Q_(ptf.Q_),
    omega_(ptf.omega_)
{}


Foam::myQFourierFuncFvPatchVectorField::
myQFourierFuncFvPatchVectorField
(
    const myQFourierFuncFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    Q_(ptf.Q_),
    omega_(ptf.omega_)
{}


Foam::myQFourierFuncFvPatchVectorField::
myQFourierFuncFvPatchVectorField
(
    const myQFourierFuncFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    Q_(ptf.Q_),
    omega_(ptf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::myQFourierFuncFvPatchVectorField::autoMap
// (
//     const fvPatchFieldMapper& m
// )
// {
//     fixedValueFvPatchVectorField::autoMap(m);
//     fieldData_.autoMap(m);
// }
//
//
// void Foam::myQFourierFuncFvPatchVectorField::rmap
// (
//     const fvPatchVectorField& ptf,
//     const labelList& addr
// )
// {
//     fixedValueFvPatchVectorField::rmap(ptf, addr);
//
//     const myQFourierFuncFvPatchVectorField& tiptf =
//         refCast<const myQFourierFuncFvPatchVectorField>(ptf);
//
//     fieldData_.rmap(tiptf.fieldData_, addr);
// }


void Foam::myQFourierFuncFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get patch range and orientation
    boundBox bb(patch().patch().localPoints(), true); // true stands for (doReduce = true) retaled with parallel case

    vector ctr = 0.5*(bb.max() + bb.min()); // center point of the patch

    const vectorField& c( (patch().Cf()) ); // coordinate of each face center on the patch

    scalarField r_(mag(c - ctr)); // Calculate the radius of each face from the center

    R_ = gMax(r_); // maximum radius

    tmp<vectorField> valuesV((-1)*u_t(r_, R_, t())*patch().nf());

    fixedValueFvPatchVectorField::operator==
    (
    valuesV
    );


    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::myQFourierFuncFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("Q",Q_);
    os.writeEntry("omega",omega_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myQFourierFuncFvPatchVectorField
    );
}

// ************************************************************************* //
