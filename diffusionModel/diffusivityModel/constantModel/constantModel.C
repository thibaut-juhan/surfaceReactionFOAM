/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "constantModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Constructeur * * * * * * * * * * * //

Foam::constantModel::constantModel
(
	volScalarField& T,
	volScalarField& p,
	dictionary& dict,
	PtrList<volScalarField>& Y
):basicDiffusivityModel(T,p,dict),Y_(Y)
{
	for(int i = 0;i < Y_.size() ; i++)
	{
		Dval_[i]  = readScalar(dict.subDict("diffusionCoefficient").lookup("D"+Y_[i].name()));
	}
}
// *************************************************** //

Foam::PtrList<Foam::volScalarField> Foam::constantModel:: D() const
{
	const fvMesh& mesh = this->T_.mesh(); 
	PtrList<volScalarField> td(Y_.size());
	forAll(td,i)
	{
		td.set
		(
		 	i,
			new volScalarField
			(
				IOobject
				(
					"D" + Y_[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar(dimViscosity,Zero)
			)
   		);
	}
	
	forAll(Y_,i)
	{
		td[i]  = this->Dval_[i];
	}

	//-------------- COEFFICIENT DE DIFFUSION SUR LES BC --------//
	forAll(T_.boundaryField(),patchI)
	{	
		for(int i = 0 ; i < Y_.size() ; i++)
		{
			fvPatchScalarField& D_patch = td[i].boundaryFieldRef()[patchI];
			D_patch = this->Dval_[i];
		}
	}

	
	return td;
}
