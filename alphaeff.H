volTensorField kappaeffbycp = kappaeff / thermo.Cp();

// Overwrite alphaeff tensor to the corresponding fluid 
// thermal diffusivity for zones without catalyst
volScalarField kappa_f = thermo.kappa();
volScalarField cp_f    = thermo.Cp();

forAll(cellsNoCatalyst, cc)
{
    label celli          = cellsNoCatalyst[cc];
    scalar kappafbycpf   = kappa_f[celli] / cp_f[celli];
    kappaeffbycp[celli]  = tensor( kappafbycpf, 0, 0, 0, kappafbycpf, 0, 0, 0, kappafbycpf );
}

// Check whether the conductivity is zero for any diagonal member in kappaeff
// The zero component is the one indicating the streamwise direction
// In such direction, the effective thermal conductivity is calculated by a 
// parallel-like circuit for the solid and fluid thermal conductivities
// That is, keff = ksolid*(1-epsilon) + kfluid*epsilon
// In the other directions, the circuit is series-like
// That is, keff = 1.0/[ (1-epsilon)/ksolid + epsilon/kfluid ]
// Note that ksolid is uniform
label refcell = 0;
if( kappaeff[refcell].component(tensor::ZZ) == 0.0 )
{
	forAll(cellsCatalyst, cc)
	{
    	label celli                               = cellsCatalyst[cc];
        scalar kappa_s                            = kappaeff[refcell].component(tensor::YY);
    	kappaeffbycp[celli].component(tensor::ZZ) = ( kappa_f[celli]*epspore[celli] + kappa_s*(1.0 - epspore[celli]) )/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::XX) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::YY) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
    }
}
else if( kappaeff[refcell].component(tensor::XX) == 0.0 )
{
	forAll(cellsCatalyst, cc)
    {
        label celli                               = cellsCatalyst[cc];
        scalar kappa_s                            = kappaeff[refcell].component(tensor::YY);
        kappaeffbycp[celli].component(tensor::ZZ) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::XX) = ( kappa_f[celli]*epspore[celli] + kappa_s*(1.0 - epspore[celli]) )/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::YY) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
    }
}
else if( kappaeff[0].component(tensor::YY) == 0.0 )
{
    forAll(cellsCatalyst, cc)
    {
        label celli                               = cellsCatalyst[cc];
        scalar kappa_s                            = kappaeff[refcell].component(tensor::XX);
        kappaeffbycp[celli].component(tensor::ZZ) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::XX) = (1.0/( epspore[celli]/kappa_f[celli] + (1.0 - epspore[celli])/kappa_s ))/ cp_f[celli];
        kappaeffbycp[celli].component(tensor::YY) = ( kappa_f[celli]*epspore[celli] + kappa_s*(1.0 - epspore[celli]) )/ cp_f[celli];
    }
}