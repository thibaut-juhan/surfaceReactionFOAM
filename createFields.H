#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");
basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y      = composition.Y();
const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();
#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// create a reaction model -> reaction en volume 
Info<< "Creating reaction model\n" << endl;
autoPtr<BasicChemistryModel<psiReactionThermo>> pChemistry
(
    BasicChemistryModel<psiReactionThermo>::New(thermo)
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());
// *************** UN SEUL COEFFICIENT DE DIFFUSION ****************** //
IOdictionary speciesDict
(   
	IOobject
	(
		"speciesDiffusionProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    	)
);
wordList speciesName
(
    speciesDict.lookup("species")
);

//-*******************  DIFFUSION COEFFICIENT *******************-//

Info << "Create list of specie molar weight.." << endl;
List<scalar> Wi (Y.size());
for (int i = 0; i < Y.size();i++)
{
	Wi[i] = readScalar(speciesDict.subDict("molarWeight").lookup("M"+Y[i].name()));
}	

Info << "Create fileds for diffusion coefficient.."<<endl;
PtrList<volScalarField> Di(Y.size());
forAll(Di,i)
{
	Di.set(i,new volScalarField
	              (	
		       	 IOobject
			 (
				"D" +  speciesName[i],
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			 ),
			 mesh,
    			 dimensionedScalar(dimViscosity,Zero)	 
		      )
       	       );
}

if(word(speciesDict.lookup("method")) == "average")
{
	Info << "average mixture model .. " << endl;
}
else if (word(speciesDict.lookup("method")) == "constant")
{
	Info << "constant model .. " << endl;
}

//-******************************************************************-//
/*
Info<< "Create field for viscosity.." << endl;
volScalarField mu
(
	IOobject
	(
		"muEff",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar(dimDynamicViscosity,Zero)
);
*/
//-********************** SPECIFIC ENTHALPY ***********************************-//
Info<<"Create fields for specific enthalpy.. \n"<< endl;
PtrList<volScalarField> h_i (Y.size());
forAll(h_i,i)
{
	h_i.set(i,new volScalarField
	              (	
		       	 IOobject
			 (
				"hi",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			 ),
			 mesh,
    			 dimensionedScalar(dimEnergy/dimMass,Zero)	 
		      )
       	       );
}
//taux de reaction de l'espece i 
Info<<"Create fields for reaction rate ..\n"<< endl;
PtrList<volScalarField> RR_i (Y.size());
forAll(RR_i,i)
{
	RR_i.set(i,new volScalarField
	              (	
		       	 IOobject
			 (
				"RR_" + speciesName[i],
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			 ),
			 mesh,
    			 dimensionedScalar(dimDensity/dimTime,Zero)	 
		      )
       	       );
}
Info<<"Create diffusive flux.. \n"<<endl;
PtrList<surfaceScalarField> Ji(Y.size()); // flux defini en tant que surfaceScalarField
forAll(Ji,i)
{
	Ji.set(i,new surfaceScalarField
	              (	
		       	 IOobject
			 (
				"J"+speciesName[i],
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			 ),
			 mesh,
    			 dimensionedScalar(dimMass/dimTime,Zero)	 
		      )
       	       );
}
// flux correctif -> la vitesse de diffusion Vi_corrige = Vi + Vc -> modifie le flux
surfaceScalarField Jc
(
 	IOobject
	(
		"Jc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar(dimMass/dimTime,Zero)
);
// *********************************************************************** //
volScalarField Qdot
(   
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);
// ******************************************************************* //
#include "createDpdt.H"
#include "createK.H"
#include "createMRF.H"
#include "createFvOptions.H"
