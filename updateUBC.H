
//-------------------- MODIFIE LA CONDITION AU LIMITE SUR CATALYTIC WALL --------------//
const label patchID               = mesh.boundaryMesh().findPatchID("catalyticWall"); 
const fvPatchScalarField rhoPatch = rho.boundaryField()[patchID];
fvPatchVectorField& Unormal       = U.boundaryFieldRef()[patchID];
scalarField sumRRi(T.boundaryField()[patchID].size());
vectorField normalFace = mesh.Sf().boundaryField()[patchID]/(mesh.magSf().boundaryField()[patchID]);
if(U.boundaryField()[patchID].type() == fixedValueFvPatchVectorField::typeName)
{
	forAll(rhoPatch,faceI)
	{
		/*
		const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
		const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
		*/

		for(int i = 0;i<Y.size();i++)
		{
			fvPatchScalarField& RRPatchI = RR_i[i].boundaryFieldRef()[patchID];
			sumRRi[faceI]               += RRPatchI[faceI];	
		}
		// parcours le vecteur :
		for(int j = 0;j<normalFace;j++)
		{
			if(normalFace[faceI][j] != 0.)
			{
				Unormal[faceI][j] = sumRRi[facei]/rhoPatch[faceI];
			}
			else
			{
				Unormal[faceI][j] = 0.0;
			}
		}
		
	}	
}


