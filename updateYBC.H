// #include "updateUBC.H" -> vitesse de Stefan : hypotheses non-glissement
// En prenant en compte la vitesse de Stefan -> necessite une correction sur
// la condition aux limites sur Yi
{
    //-------------------- MODIFIE LA CONDITION AU LIMITE SUR CATALYTIC WALL --------------//
    const label patchID               = mesh.boundaryMesh().findPatchID("catalyticWall");
    const fvPatchScalarField rhoPatch = rho.boundaryField()[patchID]; 
    scalarField Yt(0.0*Y[0].boundaryField()[patchID]);
    for(int i =0; i<Y.size(); i++)
    {
	volScalarField& Yi = Y[i];
	if( Yi.boundaryField()[patchID].type() == fixedGradientFvPatchScalarField::typeName
			&& (i!= inertIndex && composition.active(i)) )
	{
		// boundaryFieldRef() -> reference non cste
		// boundaryField()    -> reference cste

		fvPatchScalarField  DiPatch  = Di[i].boundaryField()[patchID];
		fvPatchScalarField  RRiPatch = RR_i[i].boundaryField()[patchID];
	    	fixedGradientFvPatchScalarField& Ygrad = refCast<fixedGradientFvPatchField<scalar> >(
				Yi.boundaryFieldRef()[patchID]);
		scalarField& YpatchGrad = Ygrad.gradient();
		forAll(RRiPatch,faceI)
		{
			YpatchGrad[faceI] = effSurface*( RRiPatch[faceI]/(DiPatch[faceI]*rhoPatch[faceI]) );
			if(Yi.boundaryField()[patchID][faceI] < 0.)
			{
					Pout << "WARNING : Yi<0" << endl;
			}
		}
	}
    }
}	
