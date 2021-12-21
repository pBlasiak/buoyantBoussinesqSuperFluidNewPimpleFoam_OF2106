/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    GalusinskiVigneauxNum

Description
    Calculates and outputs the mean and maximum GalusinskiVigneaux Numbers.

References:
	[1] T. Suekane, M. Sekiguichi, S. Hirai, T. Okamura, 2003 
	    Heat transfer and flow of He II in narrow channels, Cryogenics, 45, 125-131


\*---------------------------------------------------------------------------*/

scalar SuekaneNum = 0.0;
//scalar meanGalusinskiVigneauxNum = 0.0;

    SuekaneNum = gMax(3*runTime.deltaTValue()/
			(pow(laminarTransport.onebyf()*magGradT2, -1./3)*laminarTransport.rhoHe()
	*laminarTransport.cpHe()*pow(CellDims.field(),2)));


    //meanGalusinskiVigneauxNum = gAverage( 1.0/(0.5*(C2*tauMu + sqrt(pow(C2*tauMu, 2) + 4*C1*tauRho*tauRho)))*
//								runTime.deltaTValue());

//Info<< "GalusinskiVigneaux Number mean: " << meanGalusinskiVigneauxNum
   Info << "SuekaneNum max: " << SuekaneNum << endl;
   //Info << "delta t powinno byc: " << gMax(pow(laminarTransport.onebyf()*magGradT2, -1./3)*laminarTransport.rhoHe()
	//*laminarTransport.cpHe()*pow(CellDims.field(),2)/3) << endl;
//Info<< "GV dt: " << GVdt << endl;

//// Courant number calculation for comparison
//scalar CoNum = 0.0;
//scalar meanCoNum = 0.0;
//
//if (mesh.nInternalFaces())
//{
//    scalarField sumPhi
//    (
//        fvc::surfaceSum(mag(phi))().primitiveField()
//    );
//
//    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
//
//    meanCoNum =
//        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
//}
//
//Info<< "Courant Number mean: " << meanCoNum
//    << " max: " << CoNum << endl;
// ************************************************************************* //
