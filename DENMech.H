if( deactivationChemistry) 
{

	// Gas constant units in J/(kmol * K)
    const scalar RuSI  = 8314.463;

	// DEN factor is dimensionless
	forAll(cellsCatalyst, cc){
		label celli = cellsCatalyst[cc];
		DENfact[celli] = 1.0;
	}

	// Auxiliary array for molecular weight
	volScalarField MW = thermo.W();

	// Auxiliary array for H2 composition
	volScalarField& YH2 = Y[indexH2ref];

	// DENConc is dimensionless
	forAll(cellsCatalyst, cc){
		label celli = cellsCatalyst[cc];
		DENConc[celli] = rho[celli]*YH2[celli]/composition.Wi(indexH2ref);
	}
	
    forAll(Y, i)
    {

        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            // Units of E0 in J/kmol
            scalar energyFormation = H0f[i];

			if(i != dimensionlessKIndex )
			{
				
				// Pre-exponential factor units Pa^(-1) [Reaction 0 and 2]
        		scalar AKi = adsorptionKFactors[i];

        		if( AKi > 0.0 ){
					forAll(cellsCatalyst, cc)
					{
						label celli = cellsCatalyst[cc];
						// NOTE: composition.Wi(i) returns a nondimensional
        				// molecular weight
        				// Also, thermo.W() returns the molecular weight in kg/kmol
        				DENfact[celli] += AKi
        					*Foam::exp(-energyFormation/(RuSI*T[celli]))
        					*p[celli]*Yi[celli]*MW[celli]/composition.Wi(i);
					}
				}
    		}
    		else
    		{
        		// Dimensionless pre-exponential factor
        		scalar AKi = adsorptionKFactors[i];

				if( AKi > 0.0 ){
        			forAll(cellsCatalyst, cc)
					{
						label celli = cellsCatalyst[cc];
						// NOTE: composition.Wi(i) returns a nondimensional
        				// molecular weight
        				// Also, thermo.W() returns the molecular weight in kg/kmol
        				DENfact[celli] += AKi
        					*Foam::exp(-energyFormation/(RuSI*T[celli]))
        					*Yi[celli]/composition.Wi(i) 
        					/ ( YH2[celli]/composition.Wi(indexH2ref) );
					}
				}
    		}

    	}
    }
	
    DENfact            = DENfact*DENfact;
    
}
else
{
	// calculate regular reaction rate
	chemistry.calculate();
}