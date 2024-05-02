{
    //-------------------- MODIFIE LA CONDITION AU LIMITE SUR CATALYTIC WALL --------------//
    volScalarField& T              = thermo.T();
    volScalarField  Cp             = thermo.Cp();
    volScalarField  alpha          = turbulence->alphaEff();

    const label patchID            = mesh.boundaryMesh().findPatchID("catalyticWall");
    const fvPatch& patchCW         = mesh.boundary()[patchID];

    fvPatchScalarField   alphaPatch = alpha.boundaryField()[patchID];
    fvPatchScalarField&  Tpatch     = T.boundaryFieldRef() [patchID];
    fvPatchScalarField   QdotPatch  = Qdot.boundaryField() [patchID];
    fvPatchScalarField   CpPatch    = Cp.boundaryField()   [patchID]; 
    fvPatchScalarField   rhoPatch   = rho.boundaryField()  [patchID];
    
    scalarField          deltaN     = 1./mesh.deltaCoeffs().boundaryField()[patchID];
    forAll(QdotPatch,faceI)
    {
	const scalar x         = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
	const scalar y         = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
	const label faceCellsI = patchCW.faceCells()[faceI];
	const scalar lambda    = alphaPatch[faceI]*CpPatch[faceI]*rhoPatch[faceI];
	// Reaction exothermique : modifie la température sur la surface active
	Tpatch[faceI] = T[faceCellsI] + deltaN[faceI] * ( QdotPatch[faceI]/lambda );
    }

}	
