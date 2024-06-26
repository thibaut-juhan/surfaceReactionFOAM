/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::basicChemistryModel

Description
    Base class for chemistry models

SourceFiles
    basicChemistryModelI.H
    basicChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicChemistryModel_H
#define basicChemistryModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class basicChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basicChemistryModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        basicChemistryModel(const basicChemistryModel&);

        //- No copy assignment
        void operator=(const basicChemistryModel&) = delete;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Chemistry activation switch
        Switch chemistryVol_;
	Switch chemistrySurf_;

	//- Catalytic sample length
	scalar catalyticLength_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField deltaTChem_;

        // Protected Member Functions

        //- Return non-const access to the latest estimation of integration
        //  step, e.g. for multi-chemistry model
        inline volScalarField& deltaTChem();

        //- Correct function - updates due to mesh changes
        void correct();

public:

    //- Runtime type information
    TypeName("basicChemistryModel");

    // Constructors

        //- Construct from thermo
        basicChemistryModel(basicThermo& thermo);


    // Selectors

        //- Generic New for each of the related chemistry model
        template<class ChemistryModel>
        static autoPtr<ChemistryModel> New
        (
            typename ChemistryModel::reactionThermo& thermo
        );


    //- Destructor
    virtual ~basicChemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Chemistry activation switch
        inline Switch chemistryVol()  const;
	inline Switch chemistrySurf() const;

	//- catalytic sample length : 
	inline scalar catalyticLength() const;

        //- The number of species
        virtual label nSpecie() const = 0;

        //- The number of reactions
        virtual label nReaction() const = 0;

        //- Return the latest estimation of integration step
        inline const volScalarField& deltaTChem()     const;

        // Functions to be derived in derived classes

            // Fields

                //- Return const access to chemical source terms [kg/m3/s]
                virtual const volScalarField& RR
                (
                    const label i
                ) const = 0;

                //- Return access to chemical source terms [kg/m3/s]
                virtual volScalarField& RR
                (
                    const label i
                ) = 0;
                //- Return reaction rate of the speciei in reactioni
                virtual tmp<volScalarField> calculateRR
                (
                    const label reactioni,
                    const label speciei
                ) const = 0;

            // Chemistry solution
                //- Calculates the reaction rates
                virtual void calculate() = 0;
                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalar deltaT) = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const volScalarField& deltaT) = 0;

                //- Return the chemical time scale
                virtual tmp<volScalarField> tc() const = 0;

                //- Return the heat release rate [kg/m/s3]
                virtual tmp<volScalarField> Qdot() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicChemistryModelI.H"

#ifdef NoRepository
    #include "basicChemistryModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
