{
    volScalarField&  he  = thermo.he();
     
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt 
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        fvOptions(rho,he)
        + Qdot
    );

    // *********************** AJOUT ************************* //	
    for ( int i = 0 ; i < h_i.size() ; i++)
    {
	for (int j = 0 ; j < h_i[i].size() ; j++)
	{
		h_i[i][j] = thermo.composition().Hs(i,p[j],T[j]);
	}
    }
    // *********** ENTHALPIE SPECIFIQUE ********************** //

    forAll(Y,i)
    // ajout de la diffusion thermique -> diffusion des especes
    {
	volScalarField& hi = h_i[i];
	EEqn   -= fvc::div(Ji[i],hi,"div(Ji,hi)"); 
    }
    // ****************************************************** //
    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
