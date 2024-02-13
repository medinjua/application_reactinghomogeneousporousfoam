if( deactivationChemistry )
{
	// First, reinitialize omegarr_i to zero
    omegarr_i *= 0.0;

    // Compare species name to necessary species label
    // Yi and i have been defined in YEqn.H
    // NOTE: Reaction rates for species are in units kg/m3/s
    // saved member ROP is, however, in kmol/m3/s
    if(Yi.member() == "H2"){

    	volScalarField::Internal RR0(chemistry.calculateRR(0, i));
    	volScalarField::Internal RR1(chemistry.calculateRR(1, i));
    	volScalarField::Internal RR2(chemistry.calculateRR(2, i));

    	forAll(cellsCatalyst, cc)
    	{
            label celli = cellsCatalyst[cc];
    		omegarr_i[celli] = ( RR0[celli]/Foam::pow(DENConc[celli], expH2DENConc[0])
    				+ RR1[celli]/Foam::pow(DENConc[celli], expH2DENConc[1])
    				+ RR2[celli]/Foam::pow(DENConc[celli], expH2DENConc[2]) )
    				/DENfact[celli];
    	}

    }
    else if(Yi.member() == "H2O"){

    	volScalarField::Internal RR0(chemistry.calculateRR(0, i));
    	volScalarField::Internal RR1(chemistry.calculateRR(1, i));
    	volScalarField::Internal RR2(chemistry.calculateRR(2, i));

    	forAll(cellsCatalyst, cc)
    	{
            label celli = cellsCatalyst[cc];
    		omegarr_i[celli] = ( RR0[celli]/Foam::pow(DENConc[celli], expH2DENConc[0])
    				+ RR1[celli]/Foam::pow(DENConc[celli], expH2DENConc[1])
    				+ RR2[celli]/Foam::pow(DENConc[celli], expH2DENConc[2]) )
    				/DENfact[celli];
    	}

    }
    else if(Yi.member() == "CH4"){

    	volScalarField::Internal RR0(chemistry.calculateRR(0, i));
    	volScalarField::Internal RR2(chemistry.calculateRR(2, i));

    	forAll(cellsCatalyst, cc)
    	{
            label celli = cellsCatalyst[cc];
    		omegarr_i[celli] = ( RR0[celli]/Foam::pow(DENConc[celli], expH2DENConc[0])
    				+ RR2[celli]/Foam::pow(DENConc[celli], expH2DENConc[2]) )
    				/DENfact[celli];
    	}

    }
    else if(Yi.member() == "CO"){

    	volScalarField::Internal RR0(chemistry.calculateRR(0, i));
    	volScalarField::Internal RR1(chemistry.calculateRR(1, i));

    	forAll(cellsCatalyst, cc)
    	{
            label celli = cellsCatalyst[cc];
    		omegarr_i[celli] = ( RR0[celli]/Foam::pow(DENConc[celli], expH2DENConc[0])
    				+ RR1[celli]/Foam::pow(DENConc[celli], expH2DENConc[1]) )
    				/DENfact[celli];

            // For the first reaction, the default reference specie is CO
            // save the RR of this specie as the reference one
            ROP_[0][celli]    = RR0[celli]/(Foam::pow(DENConc[celli], expH2DENConc[0])
                *DENfact[celli]*composition.Wi(i));
    	}
        
    }
    else if(Yi.member() == "CO2"){

    	volScalarField::Internal RR1(chemistry.calculateRR(1, i));
    	volScalarField::Internal RR2(chemistry.calculateRR(2, i));

    	forAll(cellsCatalyst, cc)
    	{
            label celli = cellsCatalyst[cc];
    		omegarr_i[celli] = ( RR1[celli]/Foam::pow(DENConc[celli], expH2DENConc[1])
    				+ RR2[celli]/Foam::pow(DENConc[celli], expH2DENConc[2]) )
    				/DENfact[celli];

            // For the second and third reaction, the default reference specie is CO2
            // save the RR of this specie as the reference one
            ROP_[1][celli]    = RR1[celli]/(Foam::pow(DENConc[celli], expH2DENConc[1])
                *DENfact[celli]*composition.Wi(i));
            ROP_[2][celli]    = RR2[celli]/(Foam::pow(DENConc[celli], expH2DENConc[2])
                *DENfact[celli]*composition.Wi(i));
    	}

    }

    // Modify omegarr_i (multiply by catalyst density)
    // Currently, it is in kg/kgcat/s, it is necessary to multiply by the catalyst density
    // in order to get it as kg/m3/s
    omegarr_i  *= rhocat;

}
else{

    // regular reaction matrix
    // reaction rate obtained by chemistry.RR(i), previously set by chemistry.calculate()
    // NOTE: Reaction rates for species are in units kg/m3/s
    volScalarField::Internal RRi(chemistry.RR(i));

    // First, reinitialize omegarr_i to zero
    omegarr_i *= 0.0;
    forAll(cellsCatalyst, cc)
    {
        label celli       = cellsCatalyst[cc];
        omegarr_i[celli]  = RRi[celli];
    }

}