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

#include "mixtureAverageModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Constructeur * * * * * * * * * * * //

Foam::mixtureAverageModel::mixtureAverageModel
(
	Foam::volScalarField& T,
	Foam::volScalarField& p,
	Foam::dictionary& specieDict,
 	Foam::PtrList<Foam::volScalarField>& Y,
      	Foam::List<Foam::scalar>& W

):basicDiffusivityModel(T,p,specieDict),Y_(Y),W_(W)
{
}

Foam::volScalarField Foam::mixtureAverageModel::Wmean() const
{

	const fvMesh& mesh = T_.mesh();
	volScalarField Wmean
	(
		IOobject
		(
			"Dtemp",
			mesh.time().timeName(),
			mesh,
	 		IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
	 	dimensionedScalar(dimensionSet(0,0,0,0,0,0,0),Zero)
	);
	
	for(int i = 0 ; i <W_.size(); i++)
	{
		Wmean += Y_[i]*W_[i];
	}
	
	// condition aux limites:
	forAll(T_.boundaryField(),patchI)
	{
		fvPatchScalarField& WmeanPatch = Wmean.boundaryFieldRef()[patchI];
		for(int i = 0; i < W_.size(); i++)
		{
			fvPatchScalarField YiPatch     = Y_[i].boundaryField()[patchI];
			forAll(WmeanPatch,faceI)
			{
				WmeanPatch[faceI] += YiPatch[faceI]*W_[i];
			}
		}
	}

	return Wmean;
}

// coefficient de diffusion moyen dans le mélange
Foam::PtrList<Foam::volScalarField> Foam::mixtureAverageModel:: D () const
{
	const fvMesh& mesh = T_.mesh();
	PtrList<volScalarField> Deff (Y_.size());
	for(int i = 0;i<Y_.size();i++)
	{
		Deff.set
		(
		 	i,
	 		new volScalarField
			(
				IOobject
				(
					"Deff" + Y_[i].name(),
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

	volScalarField Dtemp
	(
		IOobject
		(
			"Dtemp",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar(dimTime/dimArea,Zero)
	);


	volScalarField Wmean = this->Wmean();

	for(int i = 0; i < Y_.size() ; i++)
	{	
		for(int j = 0; j < Y_.size()  ; j++)
		{
			if(i != j)
			{
				chapmanEnskogModel Dij
				(		 
					T_,
					p_,
					specieDict_, // nom du dict
					Y_[i].name(),
					Y_[j].name()
				);

				Dtemp += Wmean*Y_[j]/(W_[j]*Dij.D());
			}
		}
		//------------- CORRECTION POUR LES ESPECES INERTES -------------//
		Deff[i] = (1.-Y_[i])/(Dtemp + dimensionedScalar("SMALL",dimTime/dimArea,Foam::SMALL));
		forAll(Y_[i],cellI)
		{
			if(1 - Y_[i][cellI] == 0)
			{
				Deff[i][cellI] = 0.;
			}
		}
		Dtemp   = dimensionedScalar(dimTime/dimArea,Zero);
	}

	//---------------------- CONDITION AUX LIMITES --------------- //
	forAll(T_.boundaryField(),patchI)
	{
		fvPatchScalarField  WmeanPatch = Wmean.boundaryField()[patchI];
		fvPatchScalarField  DtempPatch = Dtemp.boundaryField()[patchI];
		for(int i = 0; i< Y_.size();i++)
		{
			fvPatchScalarField& DeffPatch  = Deff[i].boundaryFieldRef()[patchI];
			fvPatchScalarField  YiPatch    = Y_[i].boundaryField()[patchI];

			for(int j = 0; j < Y_.size();j++)
			{
				fvPatchScalarField YjPatch = Y_[j].boundaryField()[patchI];
				if (i !=j)
				{
					chapmanEnskogModel DijModel
					(		 
						T_,
						p_,
						specieDict_, // nom du dict
						Y_[i].name(),
						Y_[j].name()
					);

					volScalarField Dij = DijModel.D();
					fvPatchScalarField DijPatch = Dij.boundaryField()[patchI];
					DtempPatch += WmeanPatch*YjPatch/(W_[j]*DijPatch);
				}
			}

			DeffPatch = (1.-YiPatch)/( DtempPatch + Foam::SMALL );

			// correction pour l'espece inerte:
			forAll(YiPatch,faceI)
			{
				if(1 -YiPatch[faceI] == 0)
				{
					DeffPatch = 0.;
				}
			}

			DtempPatch   = Zero;
		}
	}
	
	return Deff;
}

Foam::tmp<Foam::volScalarField> Foam::mixtureAverageModel::phiAB
(
	Foam::volScalarField muA,
	Foam::volScalarField muB,
	Foam::scalar WA,
	Foam::scalar WB
) const
{

	const fvMesh& mesh = T_.mesh();
	tmp<volScalarField> tPhi
	(
		new volScalarField
		(
		 	IOobject
			(
				"phi_ij",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar( dimensionSet(0,0,0,0,0,0,0), Zero)
		)
	);

	volScalarField& phi = tPhi.ref();
	phi = pow( (1 + sqrt(muA/muB)* sqrt(sqrt(WB/WA)) ),2 )/sqrt(8*(1+WA/WB));

	//--- Conditions aux limites ---- //
	forAll(T_.boundaryField(),patchI)
	{
		fvPatchScalarField  muA_patch  = muA.boundaryField()[patchI];
		fvPatchScalarField  muB_patch  = muB.boundaryField()[patchI];
		fvPatchScalarField& phi_patch  = phi.boundaryFieldRef()[patchI];
			
		phi_patch = pow( (1 + sqrt(muA_patch/muB_patch)*
					sqrt(sqrt(WB/WA)) ),2 )/sqrt(8*(1+WA/WB));
	}

	return tPhi;
}

Foam::tmp<Foam::volScalarField> Foam::mixtureAverageModel::mu() const
{
	const fvMesh& mesh = T_.mesh();
	tmp<volScalarField> tMu
	(
		new volScalarField
		(
		 	IOobject
			(
				"muEff",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar(dimDynamicViscosity,Zero)
		)
	);

	volScalarField& muEff = tMu.ref();
	
	volScalarField muA
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
	);

	volScalarField muB
	(
		IOobject
		(
			"mu_j",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar(dimDynamicViscosity,Zero)
	);
	
	//- calcule la somme de Phi_ij*Yj/Wj sur j !!!
	volScalarField vTemp
	(
		IOobject
		(
			"vTemp",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar( dimensionSet(0,0,0,0,0,0,0), Zero)
	);

	for(int i = 0; i < Y_.size(); i++)
	{
		chapmanEnskogModel viscosityModelI
		(
			T_,
			p_,
			specieDict_, // nom du dict
			Y_[i].name(),
			Y_[i].name()
		);
	
		muA = viscosityModelI.mu();

		for(int j = 0; j< Y_.size(); j++)
		{
			chapmanEnskogModel viscosityModelJ
			(
				T_,
				p_,
				specieDict_,
				Y_[j].name(),
				Y_[i].name()
			);

			muB     = viscosityModelJ.mu(); 
			vTemp += Y_[j]/W_[j] * this->phiAB(muA,muB,W_[i],W_[j]);
		}

		muEff += (Y_[i]/W_[i] * muA)/vTemp;
	}
	
	return tMu;
}






//-******************************************************************-//
