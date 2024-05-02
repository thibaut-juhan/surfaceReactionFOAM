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
   Coefficient de diffusion constant -> indépendant de la température + saisi par l'utilisateur
SourceFiles
    constantModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantModel_H
#define constantModel_H

#include "basicDiffusivityModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class chapmanEnskogModel Declaration
\*---------------------------------------------------------------------------*/

class  constantModel : public basicDiffusivityModel
{
protected: 
	Foam::PtrList<Foam::volScalarField>& Y_;
	Foam::List<scalar> Dval_; // coefficient de diffusion -> lu dans les fichiers
public:
	
	// constructeur
	constantModel
	(
		Foam::volScalarField& T,
		Foam::volScalarField& p,
		Foam::dictionary& specieDict,
		Foam::PtrList<Foam::volScalarField>& Y
	);

	//-destructeur
	virtual ~constantModel() = default;
	//-fonction -> retourne le coefficient de diffusion D:
	Foam::PtrList<Foam::volScalarField> D() const; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
