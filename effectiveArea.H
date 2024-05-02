// ------------------------- CALCUL LA SURFACE ACTIVE EFFECTIVE ----------------------- //
// -------------------------  Aeff = Aactive/Atotale            ----------------------- //  	
scalar surfaceCatWall     = 0.; // surface totale de la paroi catalytic
scalar surfaceAct         = 0.; // surface active 
scalar effSurface         = 0.;
const label patchID       = mesh.boundaryMesh().findPatchID("catalyticWall");
forAll(mesh.boundaryMesh()[patchID],faceI)
{
	const scalar x  = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
	const scalar y  = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
	surfaceCatWall += mesh.magSf().boundaryField()[patchID][faceI];
	if(mag(x) <= chemistry.catalyticLength() && mag(y) <= chemistry.catalyticLength())
	{
		surfaceAct += mesh.magSf().boundaryField()[patchID][faceI];
	}
}
// parallelisé -> faire la somme entre les procs :
reduce(surfaceAct,sumOp<scalar>());
reduce(surfaceCatWall,sumOp<scalar>());
effSurface = surfaceAct/surfaceCatWall;

