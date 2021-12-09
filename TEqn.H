{
	// TODO: trzeba sprawdzic wallFunction alphatJayatille....
	// czy dobrze jest ten zmienny Pr
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

	alphaEff = kHe/(rho/rhok)/cp + alphat;
	//Info<< "alphaEff w TEqn.H = " << alphaEff << endl;
	//Info<< "alphat w TEqn.H = " << alphat << endl;
	//Info<< "rho w TEqn.H = " << rho << endl;
	//Info<< "rhok w TEqn.H = " << rhok << endl;
	//Info<< "cp w TEqn.H = " << cp << endl;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
	  + AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp
	//  + (
	//		mutualFrictionGM
	//	  ? AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp
	//      : AGM*rhon*rhos*pow(sHe/AGM/rhon*magGradT,4./3)/(rho/rhok)/cp*0
	//	)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);
    //T.correctBoundaryConditions();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;

	// beta musi byc updated
    //rhok = 1.0 - beta*(T - TRef);
}
