/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::ArrheniusReactionRate

Description
    Arrhenius reaction rate given by:

        k = A * T^beta * exp(-Ea/RT)

SourceFiles
    ArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusReactionRateEnergy_H
#define ArrheniusReactionRateEnergy_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class ArrheniusReactionRateEnergy;
Ostream& operator<<(Ostream&, const ArrheniusReactionRateEnergy&);


/*---------------------------------------------------------------------------*\
                    Class ArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusReactionRateEnergy
{
    // Private data

        scalar A_;    // facteur pré-exponentiel
        scalar beta_; // exposant T^beta
        scalar Ea_;   // ENERGIE D'ACTIVATION -> J.mol⁻1
	scalar R_;    // constante des gazs parfait (8.314 J.mol⁻1.K⁻1)
public:

    // Constructors

        //- Construct from components
        inline ArrheniusReactionRateEnergy
        (
            const scalar A,
            const scalar beta,
            const scalar Ea,
	    const scalar R
        );

        //- Construct from dictionary
        inline ArrheniusReactionRateEnergy
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "ArrheniusEnergy";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const ArrheniusReactionRateEnergy&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ArrheniusReactionRateEnergyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
