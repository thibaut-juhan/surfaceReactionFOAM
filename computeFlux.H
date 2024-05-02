// calcul le flux diffusif -> D*rhoI*grad(Yi) . n dS
for(int i = 0;i<Ji.size();i++) 
{
	volScalarField Yi = Y[i];
	Ji[i] = -(fvc::interpolate(Di[i]*rho))*(fvc::interpolate(fvc::grad(Yi)) & mesh.Sf());
}

// calcul le flux corrigé:
Jc *= 0.;
for(int i = 0;i<Ji.size();i++)
{
	Jc -= Ji[i]; // somme du flux diffusif sur TOUTES LES ESPECES : Vc
}

// modifie alors le flux calcule de CHAQUE espece -> J'i = Jc + Ji
for(int i = 0;i<Ji.size();i++)
{
	Ji[i] += Jc*fvc::interpolate(Y[i]); // sommes sur les espece
}

if(runTime.write())
{
	Jc.write();
	for(int i = 0;i<Ji.size();i++)
	{
		Ji[i].write();
	}
}

