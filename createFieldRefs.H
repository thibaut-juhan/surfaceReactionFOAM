const volScalarField& psi = thermo.psi();
const volScalarField& T   = thermo.T();
BasicChemistryModel<psiReactionThermo>& chemistry = pChemistry(); 
const label inertIndex(composition.species()[inertSpecie]);
