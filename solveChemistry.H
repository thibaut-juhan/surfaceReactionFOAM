// calcul les termes sources : 
{
	
	/*
	for(int i = 0;i<Y.size();i++)
	{
		singleConstantModel diffModel
		(
			thermo.T(),
			thermo.p(),
			speciesDict,
			speciesName[i]
		);
		D_i[i].dimensions().reset(diffModel.getD().dimensions());
		D_i[i] = diffModel.getD();
	}
	*/

	//solve chemistry :
	chemTime = t0;
	Info<<"Solve chemistry..."<<endl;
	while(tf>chemTime)
	{
		//Correction : si Yi<0 sur le patch -> Fixe à 0
		scalar dTChem   = chemistry.solve(1e-7);
		// Stocke le taux de reaction de chaque espece:
		for(int i = 0;i<Y.size();i++)
		{
			if(i != inertIndex && composition.active(i))
			{
				RR_i[i] = chemistry.RR(i);
			}	
		}
		Qdot     = chemistry.Qdot();
		chemTime = chemTime + dTChem;
	}
/*	
	const fvPatch& catalyticWallPatch = mesh.boundary()[patchID];
	for(int i = 0;i<Y.size();i++)
	{
		volScalarField&    Yi       =  Y[2];
		fvPatchScalarField YiPatch  =  Y[2].boundaryField()[patchID];
		fvPatchScalarField RRiPatch =  RR_i[2].boundaryField()[patchID];
		forAll(YiPatch,faceI)
		{
			const scalar x        = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y        = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			const label faceCellI = catalyticWallPatch.faceCells()[faceI];
			if( mag(x)<=Lsample.value() && mag(y) <= Lsample.value() && i != inertIndex)
			{
				Pout << "Y"      << speciesName[2] << " : " << Yi[faceCellI]   << endl;
				Pout << "Ypatch" << speciesName[2] << " : " << YiPatch[faceI]  << endl;
				Pout << "RRi"    << speciesName[2] << " : " << RRiPatch[faceI] << endl;
			}
		}

	}
*/	
	if(runTime.write())
	{
		for(int i = 0;i<Y.size();i++)
		{
			RR_i[i].write();
		}
	}
}
