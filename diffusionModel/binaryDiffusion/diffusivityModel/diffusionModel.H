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
    Foam::diffModel

Description
    abstract class defining the concept of binary diffusion coefficient

SourceFiles
    binaryModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionModel_H
#define diffusionModel_H

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

class diffusionModel{
protected:
	
	// Thermodynamics variable
	Foam::volScalarField& T_;           // Temperature
	Foam::volScalarField& p_;           // presssure
	Foam::dictionary&     dict_;        // dictionnary 
	const Foam::word&           namei_;       // nom de l'espece i
	const Foam::word&           namej_;       // nom de l'espece j
public:
	
	// constructeur
	diffusionModel(Foam::volScalarField& T, Foam::volScalarField& p, 
			Foam::dictionary& dict, 
			const Foam::word& namei, const Foam::word& namej);

	// destructeur
	virtual ~diffusionModel() = default;
	// methode virtuelle pure -> renvoie le coefficient de diffusion 
	virtual Foam::tmp<Foam::volScalarField> D   () const = 0; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
