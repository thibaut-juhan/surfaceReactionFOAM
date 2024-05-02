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
    Foam::chapmanEnskogModel

Description : 
    Chapman-Enskog theory for binary diffusion coefficient:
    Dij = 0.00266 T**3/2 / (p M_ij**0.5 sigma_ij**2 omega_D)
    -> Sigma_D : integral des collisions des particules 
    calculée à partir du Developpement de Neufeld(1972):
    A/(T')**B + C/exp(FT') + E/exp(FT') + G/exp(HT')
    avec T' = kT/epsi_ij -> epsi_ij : energie caracteristiques de Leonard-Jones
SourceFiles
    chapmanEnskogModel.C

\*---------------------------------------------------------------------------*/

#ifndef chapmanEnskogModel_H
#define chapmanEnskogModel_H

#include "diffusionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class chapmanEnskogModel Declaration
\*---------------------------------------------------------------------------*/

class chapmanEnskogModel : 
public diffusionModel
{
protected: // protected -> heritage
	
	// variables liees au modele
	Foam::scalar           a1,a2,a3,a4;  // parametres pour le potentiel Lennard-Jones
	Foam::scalar           b1,b2,b3,b4;  // parametres pour le potentiel Leannard-Jones

	Foam::scalar           W_ij;              // masse molaire de l'espece i
	Foam::scalar           sigma_ij;          // longueur caractéristique j (A°)
	Foam::scalar           epsi_ij;

public:
	
	// constructeur
	chapmanEnskogModel
	(
		Foam::volScalarField& T,
		Foam::volScalarField& p,
		Foam::dictionary&     dict,
		const Foam::word&           namei,
		const Foam::word&           namej
	);

	// destructeur
	virtual ~chapmanEnskogModel() = default;

	// Calcul le coefficient de diffusion binaire + viscosité d'une espece
	//
	Foam::tmp<Foam::volScalarField> D () const; 
	Foam::tmp<Foam::volScalarField> mu() const;

	Foam::scalar Tred    (const scalar T) const;
	Foam::scalar Tred2   (const scalar T) const;

	Foam::scalar omegaD  (const scalar T) const;
	Foam::scalar omegaD2 (const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
