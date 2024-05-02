{	
	//----------------------- CALCUL LE COEFFICIENT DE DIFFUSION --------------------- //

	// En fonction du modele choisi -> coefficient de diffusion différent
	if(word(speciesDict.lookup("method")) == "average")
	{
		mixtureAverageModel averageDiffusionCoeff
		(
		 	thermo.T(),
			p,
			speciesDict,
			Y,
			Wi
		);
		Di    = averageDiffusionCoeff.D(); // coefficient de diffusion 
		// mu    = averageDiffusionCoeff.mu();
	}
	else if(word(speciesDict.lookup("method")) == "constant")
	{
		constantModel constantDiffusionCoeff
		(
		 	thermo.T(),
			p,
			speciesDict,
			Y
		);
		Di = constantDiffusionCoeff.D(); // coefficient de diffusion 
	}
	
	if(runTime.write())
	{
		forAll(Di,i)
		{
			Di[i].write();
		}
	}
}
