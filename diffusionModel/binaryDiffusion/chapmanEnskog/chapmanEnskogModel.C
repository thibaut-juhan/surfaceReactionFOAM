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

#include "chapmanEnskogModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Constructeur * * * * * * * * * * * //

Foam::chapmanEnskogModel::chapmanEnskogModel
(
	volScalarField& T,
	volScalarField& p,
	dictionary&     dict,
	const word&           namei,
	const word&           namej
):diffusionModel(T,p,dict,namei,namej)
{
	//- coefficient pour calculer omegaD1
	a1  = 1.0548;
	a2  = 0.15504;
	a3  = 0.55909;
	a4  = 2.1705;

	//- coefficient pour calculer omegaD2
	b1  = 1.0413;
	b2  = 0.11930;
	b3  = 0.43628;
	b4  = 1.6041;

	// lire les données dans :
	// Dict -> espece (i ou j) -> sous-dictionnaire
	const scalar& W_i = readScalar(dict_.subDict("molarWeight").lookup("M" + namei));
	const scalar& W_j = readScalar(dict_.subDict("molarWeight").lookup("M" + namej));

	const scalar& epsi_i = readScalar(dict_.subDict("LJenergy").lookup("E" + namei));
	const scalar& epsi_j = readScalar(dict_.subDict("LJenergy").lookup("E" + namej));

	const scalar& sigma_i = readScalar(dict_.subDict("sigma").lookup("sigma" + namei));
	const scalar& sigma_j = readScalar(dict_.subDict("sigma").lookup("sigma" + namej));
	
	// calcul les coefficients à partir des données 
	sigma_ij = 0.5*(sigma_i+sigma_j);
	epsi_ij  = Foam::sqrt(epsi_i*epsi_j);
	W_ij     = 1e-3*0.5*( W_i*W_j/(W_i+W_j) )/6.022e23;	
}
// *************************************************** //


// Température réduite DE 2 ESPECES
Foam::scalar Foam::chapmanEnskogModel:: Tred(const scalar T) const// temperature réduite 
{
	return T/this->epsi_ij;
}

// Température réduite D'UNE ESPECE
Foam::scalar Foam::chapmanEnskogModel:: Tred2(const scalar T) const// temperature réduite 
{
	const scalar& epsi_i = readScalar(dict_.subDict("LJenergy").
			lookup("E" + this->namei_));
	return T/epsi_i;
}

// integral de collisions omegaD1
Foam::scalar Foam::chapmanEnskogModel:: omegaD(const scalar T) const
{
	scalar Tr     = this->Tred(T);
	scalar omegaD = pow(a1*Tr,-a2) + pow(Tr + a3,-a4);
	return omegaD;
}
// integral de collision omegaD2
Foam::scalar Foam::chapmanEnskogModel:: omegaD2(const scalar T) const // intégral de collision
{
	scalar Tr      = this->Tred(T);
	scalar omegaD2 = pow(b1*Tr,-b2) + pow(Tr + b3,-b4);
	return omegaD2;
}


// coefficient de diffusion pour le domaine -> defini au point P
Foam::tmp<Foam::volScalarField> Foam::chapmanEnskogModel:: D() const
{
	const fvMesh& mesh = this->T_.mesh(); 
	const scalar kb    = 1.380e-23;
	tmp<volScalarField> td 
	(
	 	new volScalarField
		(
			IOobject
			(
				"D_ij",
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
		scalar omegaD = this->omegaD(T_[cellI]);
		//- UNITE INTERNATIONAL : 
		d[cellI] = 3./8. * sqrt(constant::mathematical::pi*pow(kb*T_[cellI],3)
				/this->W_ij)/(p_[cellI]*pow(sigma_ij*1e-10,2)
					*omegaD*constant::mathematical::pi);

		// d[cellI]      = 0.00266*pow(pow(T_[cellI],3)*W_ij,0.5)/
		//	(p_[cellI]*pow(sigma_ij,2)*omegaD);
	}

	//- coefficient de diffusion sur les patch
	forAll(T_.boundaryField(),patchI)
	{	
		const fvPatchScalarField& T_patch = this->T_.boundaryField()[patchI];
		const fvPatchScalarField& P_patch = this->p_.boundaryField()[patchI];
		fvPatchScalarField&       D_patch = d.boundaryFieldRef()[patchI];
		forAll(T_patch,faceI)
		{
			scalar omegaD  = this->omegaD(T_patch[faceI]);

			D_patch[faceI] = 3./8. * sqrt(constant::mathematical::pi
				*pow(kb*T_patch[faceI],3)/this->W_ij)/(
				P_patch[faceI]*pow(sigma_ij*1e-10,2)*omegaD 
				*constant::mathematical::pi );

			// D_patch[faceI] = 0.00266*pow(pow(T_patch[faceI],3)*W_ij,0.5)/
			//		(P_patch[faceI]*pow(sigma_ij,2)*omegaD);
		}
	}
	
	return td;
}

Foam::tmp<Foam::volScalarField> Foam::chapmanEnskogModel::mu() const
{
	const scalar Na   = 6.022e23;
	const scalar kb   = 1.380e-23; // constante de Boltzmann
	scalar Wi  = readScalar(dict_.subDict("molarWeight").
			lookup("M"+this->namei_));
	scalar sigmai = readScalar(dict_.subDict("sigma").
			lookup("sigma" + this->namei_));
	
	const fvMesh& mesh    = this->T_.mesh(); 
	
	Wi     = Wi*1e-3;
	sigmai = sigmai*1e-10;

	tmp<volScalarField> tMu 
	(
	 	new volScalarField
		(
			IOobject
			(
				"mu_i",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar(dimDynamicViscosity,Zero)
		)
   	);
	
	volScalarField& mu =  tMu.ref();

	// calcul de mu
	forAll(T_,cellI)
	{
		//- UI
		scalar omegaD2 = this->omegaD2(T_[cellI]);
		mu[cellI] = 5./16. * sqrt( constant::mathematical::pi*T_[cellI]
				*kb*Wi/Na )/(constant::mathematical::pi
					*pow(sigmai,2)*omegaD2 );
	}

	//- coefficient de diffusion sur les patch
	forAll(T_.boundaryField(),patchI)
	{	
		const fvPatchScalarField& T_patch  = this->T_.boundaryField()[patchI];
		fvPatchScalarField&       mu_patch = mu.boundaryFieldRef()[patchI];
		forAll(T_patch,faceI)
		{	
			scalar omegaD2  = this->omegaD2(T_patch[faceI]);
			mu_patch[faceI] = 5./16. * sqrt(constant::mathematical::pi
					*T_patch[faceI]*kb*Wi/Na)/(constant::mathematical::pi
						*pow(sigmai,2)*omegaD2);
		}
	}

	return tMu;
}













//-******************************************************************-//
