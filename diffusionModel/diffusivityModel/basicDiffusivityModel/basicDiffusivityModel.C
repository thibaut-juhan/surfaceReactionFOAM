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

#include "basicDiffusivityModel.H"

// constructeur
Foam::basicDiffusivityModel::basicDiffusivityModel
(
		Foam::volScalarField& T,
		Foam::volScalarField& p,
		Foam::dictionary&     specieDict
):T_(T),p_(p),specieDict_(specieDict)
{
	nameMethod_ = word(specieDict_.lookup("method"));
}

// *************************************************** //

// coefficient de diffusion pour le domaine
Foam::PtrList<Foam::volScalarField> Foam::basicDiffusivityModel:: D() const
{
	PtrList<Foam::volScalarField> nullPtr;
	return nullPtr;
}

