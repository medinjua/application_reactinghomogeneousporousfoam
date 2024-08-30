{

    // NOTE concerning heat release rate
    // The heat release rate that appears in the total energy equation is
    // the enthalpy of reaction at standard conditions (deltaH_R^0/dt)
    // In the temperature equation, the heat release rate is the sum of the 
    // enthalpy of reaction with the summation of the products of species 
    // enthalpies and species reaction rates (- deltaH_R^0/dt - sum[rri*hi]). 
    // In such a case, the species enthalpies are defined, e.g., 
    // by NASA polynomials evaluated between standard temperature and the 
    // current target temperature
    // It is noted that at standard conditions, deltaH_R^0/dt is 
    // identically equal to sum[rri*hi^0], where hi^0 is the enthalpy 
    // at standard conditions
    // Therefore, for simplicity, one generally simply uses
    // HRR = - sum[rri*hi] with an h that is defined with reference 0K
    // Here, it is important to recall that the equation being solved 
    // is that for the total energy equation, using the enthalpy that is
    // correctly defined as the difference between target temperature and 
    // the standard conditions

    // Reinitialize Qdot to zero
    Qdot *= 0.0;

    if( deactivationChemistry ){
        for( label i=0; i<chemistry.nReaction(); i++ )
        {
            forAll(cellsCatalyst, cc)
            {
                label celli = cellsCatalyst[cc];
                // Get the heat release rate as -deltaH_f0/dt = sum_r ( omegarr_r * H_{f0,r} )
                // with omegarr_r in kmol/m3/s and H_f in J/kmol
                Qdot[celli] -= ROP_[i][celli]*rhocat[celli]*H0R[i];
            }
        }
    }
    else{
        
        // First, calculate the heat release rate as usual
        // Note: a call to chemistry.calculate() must have been done before
        Qdot = chemistry.Qdot().ref()*epspore;

        // Now set the HRR in the non-catalyst cells to zero
        forAll(cellsNoCatalyst, cc)
        {
            label celli = cellsNoCatalyst[cc];
            Qdot[celli]  = 0.0;
        }
    }

}
