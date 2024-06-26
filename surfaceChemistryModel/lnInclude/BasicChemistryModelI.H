/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::basicChemistryModel::mesh() const
{
    return mesh_;
}

inline Foam::Switch Foam::basicChemistryModel::chemistryVol() const
{
    return chemistryVol_;
}

inline Foam::scalar Foam::basicChemistryModel::catalyticLength() const
{
    return catalyticLength_;
}

inline Foam::Switch Foam::basicChemistryModel::chemistrySurf() const
{
    return chemistrySurf_;
}

inline const Foam::volScalarField&
Foam::basicChemistryModel::deltaTChem() const
{
    return deltaTChem_;
}

inline Foam::volScalarField&
Foam::basicChemistryModel::deltaTChem()
{
    return deltaTChem_;
}

// ************************************************************************* //
