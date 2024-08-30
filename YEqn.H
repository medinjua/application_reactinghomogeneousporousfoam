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

// Modification catalyst mechanism
#include "DENMech.H"
// End modification

// Following is commented out for explicit calculation of RR in catalyst sims
//reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        //---------------------------------------------------------------------
        // Modifications for porous catalyst simulations
        // In general, any explicit term can be added to the matrix equations
        // required is simply to insert a volScalarField (or volVectorField)
        // on the RHS of the equation (consistent units are naturally required)
        // In the following, the reaction rate term (explicit) is calculated
        // for species Yi
        #include "DENRRMech.H"     
        //---------------------------------------------------------------------
        // End modifications catalyst simulations

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          // Note: default laminar transport model for species diffusion 
          // is unityLewisFourier, the species diffusion coefficient is 
          // then simply the thermal diffusion coefficient
          // note that the documentation for divj says
          // -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi)
          // this->alpha() refers to the alphaField, simply a volume fraction
          // DEff is simply the thermal diffusion coefficient
          + thermophysicalTransport->divj(Yi)
         ==
            //reaction->R(Yi)               // replaced RR matrix
            omegarr_i
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
