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
	basicDiffusivityModel

Description
    abstract class defining the concept of specie diffusion coefficient

SourceFiles
	basicDiffusivity.C
\*---------------------------------------------------------------------------*/

#ifndef basicDiffusivityModel_H
#define basicDiffusivityModel_H

#include "specie.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "tmp.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicDiffusivityModel
{

protected:	
	// Thermodynamics variable
	Foam::volScalarField& T_;          // Temperature
	Foam::volScalarField& p_;          // presssure
	Foam::dictionary& specieDict_;
	Foam::word nameMethod_; // nom de la méthode choisi par l'utilisateur

public:
	// constructeur
	basicDiffusivityModel
	(
		Foam::volScalarField& T, 
		Foam::volScalarField& p,
	        Foam::dictionary& specieDict
	);
	// destructeur
	virtual ~basicDiffusivityModel() = default;
	// methode virtuelle pure -> renvoie le coefficient de diffusion 
	virtual Foam::PtrList<Foam::volScalarField> D() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
