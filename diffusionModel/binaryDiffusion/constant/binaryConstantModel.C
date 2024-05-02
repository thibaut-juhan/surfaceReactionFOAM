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

#include "binaryConstantModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Constructeur * * * * * * * * * * * //

Foam::binaryConstantModel::binaryConstantModel
(
	volScalarField& T,
	volScalarField& p,
	dictionary&     dict,
	word&           namei,
	word&           namej
):diffusionModel(T,p,dict,namei,namej)
{
	const word name = namei + "_" + namej;
	Dval  = readScalar(dict.subDict(
				"diffusionCoefficient").lookup("D"+name));
}
// *************************************************** //

// coefficient de diffusion pour le domaine -> defini au point P
Foam::tmp<Foam::volScalarField> Foam::binaryConstantModel:: D() const
{
	const fvMesh& mesh = this->T_.mesh(); 
	// recupere le maillage -> instancie un objet volScalarField

	tmp<volScalarField> td 
	(
		new volScalarField
		(
			IOobject
			(
				"D_" + namei_  + "_" + namej_,
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar(dimViscosity,Zero)
		)
   	);
	
	volScalarField& d = td.ref();

	// Parcours les cellules et calcul le coefficient
	forAll(T_,cellI)
	{
		d[cellI]  = this->Dval;
	}

	//- coefficient de diffusion sur les patch
	forAll(T_.boundaryField(),patchI)
	{	
		fvPatchScalarField        D_patch = d.boundaryField() [patchI];
		forAll(D_patch,faceI)
		{
			D_patch[faceI] = this->Dval;
		}
	}
	
	return td;
}
