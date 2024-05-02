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

Class
    Foam::mixtureAverageModel

Description :

    mixture average diffusion coefficient based on
    Chapman-Enskog theory for binary diffusion coefficient:
    Dij = 0.00266 T**3/2 / (p M_ij**0.5 sigma_ij**2 omega_D)
    -> Omega_D : integral des collisions des particules 
    calculée à partir du Developpement de Neufeld(1972):
    A/(T')**B + C/exp(FT') + E/exp(FT') + G/exp(HT')
    avec T' = kT/epsi_ij 
    -> epsi_ij  = sqrt(epsii.epsij) potentiel de Leonard-Jones
    -> sigma_ij = diametre des particules
SourceFiles
    mixtureAverageModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureAverageModel_H
#define mixtureAverageModel_H

#include "specie.H"
#include "PtrList.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "tmp.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "chapmanEnskogModel.H"
#include "basicDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class mixtureAverageModel Declaration
\*---------------------------------------------------------------------------*/

class mixtureAverageModel :public basicDiffusivityModel
{

protected:	
	PtrList<Foam::volScalarField>& Y_;    // liste des fractions massique Yi
	List<Foam::scalar>&  W_;    // liste des masses molaires
public:
	
	// constructeur
	mixtureAverageModel
	(
		Foam::volScalarField& T,
		Foam::volScalarField& p,
		Foam::dictionary& specieDict,
		Foam::PtrList<Foam::volScalarField>& Y,
		Foam::List<Foam::scalar>& W
	);

	// destructeur
	virtual ~mixtureAverageModel() = default;
	// fonction -> retourne le coefficient de diffusion D moyen:
	Foam::PtrList<Foam::volScalarField> D() const;
	// fonction -> retourne la viscosité dynamique moyenne dans le mélange
	Foam::tmp<Foam::volScalarField> phiAB
	(
		Foam::volScalarField muA,
		Foam::volScalarField muB,
		Foam::scalar         WA,
		Foam::scalar         WB
	) const;
	Foam::tmp<Foam::volScalarField> mu() const;
       	Foam::volScalarField Wmean() const;	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
