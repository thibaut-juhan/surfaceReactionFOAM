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
    Foam::constantModel

Description
   Coefficient de diffusion constant -> indépendant de la température
SourceFiles
    constantModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryConstantModel_H
#define binaryConstantModel_H

#include "diffusionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class chapmanEnskogModel Declaration
\*---------------------------------------------------------------------------*/

class binaryConstantModel : 
public diffusionModel
{
private: // protected -> heritage
	Foam::scalar  Dval;  // coefficient de diffusion -> lu dans les fichiers
public:
	
	// constructeur
	binaryConstantModel
	(
		Foam::volScalarField& T,
		Foam::volScalarField& p,
		Foam::dictionary&     dict,
		Foam::word&           namei,
		Foam::word&           namej
	);

	//-destructeur
	virtual ~binaryConstantModel() = default;
	//-fonction -> retourne le coefficient de diffusion D:
	Foam::tmp<Foam::volScalarField> D    () const; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
