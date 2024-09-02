{
    if(!isothermalFluid){
        //**********************************************************************
        // Modification porous solver
        // Set effective thermal diffusion coefficient tensor
        #include "alphaeff.H"
        #include "HRR.H"     
        // End modification
        //**********************************************************************

        volScalarField& he = thermo.he();

        // porous solver modifies density
        fvScalarMatrix EEqn
        (
            fvm::ddt((epspore*rho*cp_f + (1.0-epspore)*rhos*cpvs)/cp_f, he) + epspore*mvConvection->fvmDiv(phi, he)
          //+ epspore*fvc::ddt(rho, K) + epspore*fvc::div(phi, K)
          + (
                he.name() == "e"
                ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
                : -dpdt
            )
        //+ thermophysicalTransport->divq(he) // replaced to allow kappaeff
        // Note: implicit treatment with fvm::laplacian
        // explicit treatment with fvc::laplacian
          - fvm::laplacian(kappaeffbycp, he)    // see definition of
                                           // thermophysicalTransport>divq
        ==
            //reaction->Qdot()                  // replaced Qdot field by HRR
            epspore*Qdot
          + fvModels.source(rho, he)
        );

        EEqn.relax();

        fvConstraints.constrain(EEqn);

        EEqn.solve();

        fvConstraints.constrain(he);

        thermo.correct();
    }

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
