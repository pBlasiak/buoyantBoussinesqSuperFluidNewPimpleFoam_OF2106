    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

	//czy ta zmienna jest widoczna w TEqn czy znika zaraz po 
	//opuszczeniu UEqn?
	//volVectorField gradT(fvc::grad(T));
	
	// Update properties for M
	rho = laminarTransport.rhoHe();
	rhon = laminarTransport.rhon();
	rhos = laminarTransport.rhos();
	AGM = laminarTransport.AGMHe();
	sHe = laminarTransport.sHe();
	beta = laminarTransport.betaHe();
	cp = laminarTransport.cpHe();
    rhok = 1.0 - beta*(T - TRef);
	if(gradP)
	{
		Info<< "Licze z gradP!!!!" << endl;
		gradT = fvc::grad(T)-fvc::grad(p*(rho/rhok))/rho/sHe;	
	}
	else
	{
		gradT = fvc::grad(T);	
	}
	magGradT = max(mag(gradT), dimensionedScalar("small", dimTemperature/dimLength, SMALL)); 
	magGradT2 = pow(magGradT,2);
	
	// (1/(f|gradT|^2))^(1/3)
	kHe = pow(laminarTransport.onebyf()/magGradT2, 1./3);

	// To jest potrzebne do alphatJayatille....
	Pr = cp*laminarTransport.etaHe()/kHe;
	//pow(laminarTransport.onebyf()/max(pow(mag(fvc::grad(T)),2), dimensionedScalar("small", dimTemperature*dimTemperature/dimLength/dimLength, SMALL)), 1./3);

	volVectorField M
	(
		fvc::div
		(
			rhon*rhos/rho*pow(sHe/AGM/rhon/magGradT2, 2./3)*gradT*gradT
		)/(rho/rhok)
	  + laminarTransport.etaHe()/(rho/rhok)*pow
	    (
			pow(rhos,3)*sHe/AGM/pow(rho,3)/rhon/magGradT2, 1./3
		)*
		(
			fvc::laplacian(gradT) + 1./3*fvc::grad(fvc::div(gradT))
		)
	);
	Info<< "max(M) = " << max(M) << endl;

	// TODO: zastanowic sie czy nie powinno byc
	// rho pod pochodna jak jest w buoyantPimpleFoam
	// ale w sumie stosujemy aproksymacje Boussinesqa wiec jest ok:
	// In flows accompanied by heat transfer, the fluid properties are normally functions of temperature. The variations may be small and yet be the cause of the fluid motion. If the density variation is not large, one may treat the density as constant in the unsteady and convection terms, and treat it as variable only in the gravitational term. This is called the Boussinesq approximation.
	// We have to know when these approximations are valid. Ferziger and Peric (2001, p.15) states that:

    //This approximation introduces errors of the order of 1% if the temperature differences are below e.g. 2℃ for water and 15℃ for air. The error may be more substantial when temperature differences are larger; the solution may even be qualitatively wrong.
	// https://caefn.com/openfoam/solvers-buoyantboussinesqpimplefoam
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
	  //- fvm::laplacian(laminarTransport.etaHe()/(rho/rhok), U)
      //- fvc::div(laminarTransport.etaHe()/(rho/rhok)*dev(fvc::grad(U)().T()))
	  //+ laminarTransport.etaHe()/(rho/rhok)*fvm::laplacian(U)
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
				  - (fvc::interpolate(M) & mesh.Sf())
            )
        );

        fvOptions.correct(U);
    }
