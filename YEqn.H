#include "computeFlux.H"
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
{
    volScalarField Yt  (0.0*Y[0]);
    for(int i = 0; i < Y.size(); i++)
    {
        if (i != inertIndex && composition.active(i))
        {    
        	volScalarField&      Yi      = Y[i];
            	fvScalarMatrix YiEqn
            	(
                	fvm::ddt(rho, Yi)
              		+ mvConvection->fvmDiv(phi, Yi)
			- fvm::laplacian(Di[i]*rho,Yi) // flux sans correction
             		==
			RR_i[i]
			- fvm::div(Jc,Yi,"div(Jc,Yi)") // AJOUT : flux avec correction
              		+ fvOptions(rho, Yi)
            	);
            	YiEqn.relax();
            	fvOptions.constrain(YiEqn);
            	YiEqn.solve(mesh.solver("Yi"));
            	fvOptions.correct(Yi);
            	Yi.max(0.0); // assure que Yi>0
            	Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
