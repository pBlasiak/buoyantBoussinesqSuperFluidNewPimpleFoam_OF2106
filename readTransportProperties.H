singlePhaseHeliumTransportModel laminarTransport(U, phi);

// updates helium thermo-physical properties
// based on temperature
laminarTransport.correct();

//// Thermal expansion coefficient [1/K]
//dimensionedScalar beta
//(
//    "beta",
//    dimless/dimTemperature,
//    laminarTransport
//);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
//dimensionedScalar Pr("Pr", dimless, laminarTransport);
// TODO: trzeba to sie zastanowic czy dobrze jest
volScalarField Pr
(
    IOobject
    (
        "Pr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.cpHe()*laminarTransport.etaHe()/
	pow(laminarTransport.onebyf()/max(pow(mag(fvc::grad(T)),2), dimensionedScalar("small", dimTemperature*dimTemperature/dimLength/dimLength, SMALL)), 1./3)
);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

bool gradP
(
	laminarTransport.lookupOrDefault("gradP", false)
);

bool mutualFrictionGM
(
	laminarTransport.lookupOrDefault("mutualFrictionGM", false)
);
