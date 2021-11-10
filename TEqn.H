{
	// TODO: trzeba sprawdzic wallFunction alphatJayatille....
	// czy dobrze jest ten zmienny Pr
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

	alphaEff = 	pow
		(
			laminarTransport.onebyf()/max
			(
				//pow(mag(fvc::grad(T)),2), 
				pow(mag(gradT),2), 
				dimensionedScalar("small", dimTemperature*dimTemperature/dimLength/dimLength, SMALL)
			), 1./3
		)/(rho/rhok)/cp
	  + alphat
	;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
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
