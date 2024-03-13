const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

// ****************************************************************************
// Additions for porous catalyst simulations
// ****************************************************************************

// chemistry model required
basicChemistryModel& chemistry = pChemistry();

// Define the ThermoType class (required for access to reactions object and methods)
// note the defined method
// typedefThermo( sutherlandTransport, sensibleEnthalpy, janafThermo, perfectGas, specie )
typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > thermoPhysics;
standardChemistryModel<thermoPhysics> &reaction_kinetics = 
       dynamic_cast<standardChemistryModel<thermoPhysics>&>(chemistry);
PtrList<Reaction<thermoPhysics>> reactions(reaction_kinetics.reactions());

const label inertIndex(composition.species()[inertSpecie]);

// Define cell set list with and without catalyst zone
const word cellSetCatalyst = "catalyst";
label zoneID = mesh.cellZones().findZoneID(cellSetCatalyst);
if (zoneID == -1)
{
    FatalErrorIn("createFieldRefs")
        << "Cannot find cellZone " << cellSetCatalyst << endl
        << "Valid cellZones are " << mesh.cellZones().names()
        << exit(FatalError);
}
const labelList& cellsCatalyst = mesh.cellZones()[zoneID];

const word cellSetEntry = "entry";
zoneID = mesh.cellZones().findZoneID(cellSetEntry);
if (zoneID == -1)
{
    FatalErrorIn("createFieldRefs")
        << "Cannot find cellZone " << cellSetEntry << endl
        << "Valid cellZones are " << mesh.cellZones().names()
        << exit(FatalError);
}
const labelList& cellsEntry = mesh.cellZones()[zoneID];

const word cellSetExit = "exit";
zoneID = mesh.cellZones().findZoneID(cellSetExit);
if (zoneID == -1)
{
    FatalErrorIn("createFieldRefs")
        << "Cannot find cellZone " << cellSetExit << endl
        << "Valid cellZones are " << mesh.cellZones().names()
        << exit(FatalError);
}
const labelList& cellsExit = mesh.cellZones()[zoneID];

labelList tmp_cellsNoCatalyst( cellsEntry.size() + cellsExit.size() );
label     cc = 0;
forAll( cellsEntry, i ){
    tmp_cellsNoCatalyst[cc] = cellsEntry[i];
    cc += 1;
}
forAll( cellsExit, i ){
    tmp_cellsNoCatalyst[cc] = cellsExit[i];
    cc += 1;
}
const labelList& cellsNoCatalyst = tmp_cellsNoCatalyst;

// Overwrite porosity (epsilon = 1) for zones without catalyst
forAll( cellsNoCatalyst, cc ){
    label celli    = cellsNoCatalyst[cc];
    epspore[celli] = 1.0;
}

// Check for option to run isothermal simulation (deactivates energy eqn)
bool isothermalFluid = readBool(thermo.properties().lookup("isothermalFlow"));
if( isothermalFluid ){
    Info<< "\nRunning isothermalFluid mode (no energy eqn solved)" << endl;
}

// Following initialization is dummy for no deactivation chemistry
List<scalar> adsorptionKFactors(chemistry.nSpecie(), 0.0);
List<scalar> H0f(chemistry.nSpecie(), 0.0);
List<scalar> H0R(chemistry.nReaction(), 0.0);
List<scalar> expH2DENConc(chemistry.nReaction(), 0.0);
label dimensionlessKIndex(0);
label indexH2ref(0);
bool deactivationChemistry( readBool(thermo.properties().lookup("deactivationChemistry")) );
PtrList<List<specieCoeffs>> products_rrn(chemistry.nReaction());
PtrList<List<specieCoeffs>> reactants_rrn(chemistry.nReaction());

PtrList<volScalarField> ROP_(chemistry.nReaction());
for( label i=0; i<chemistry.nReaction(); i++ )
{
    scalar reactionk(i);
    const word RRname = "ROP.R" + name(reactionk);
    ROP_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                RRname,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                deactivationChemistry ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimVolume/dimTime, 0.0)
        )
    );    
}

if( deactivationChemistry ){

    Info<< "\nReading additional adsorption constant reaction data " 
        << "(deactivation chemistry, adsorptionRateConstant dictionary)" << endl;


    // Adsorption constant parameters
    IOobject adsorptionConstant
    (
        "adsorptionRateConstant",
        runTime.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ); 
    IOdictionary adsorptionK(adsorptionConstant);

    adsorptionKFactors = List<scalar>(adsorptionK.lookup("factors"));
    H0f = List<scalar>(adsorptionK.lookup("energyFormation"));
    expH2DENConc = List<scalar>(adsorptionK.lookup("expH2DENConc"));
    if (!adsorptionK.lookup("indexH2ref"))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo.properties())
            << "Reference species (H2) not found in available species "
            << "or its index is not defined in the dictionary  "
            << composition.species()
            << exit(FatalIOError);
    }
    if (!adsorptionK.lookup("dimensionlessKIndex"))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo.properties())
            << "Reference consumption species (H2O) not found in available species"
            << "or its index is not defined in the dictionary  "
            << composition.species()
            << exit(FatalIOError);
    }
    dimensionlessKIndex = readLabel(adsorptionK.lookup("dimensionlessKIndex"));
    indexH2ref = readLabel(adsorptionK.lookup("indexH2ref"));

    Info<< "\nIndex H2 is: " << indexH2ref << "  - Index reference consumption: "
        << dimensionlessKIndex << endl;

    if( gMax( Y[indexH2ref].internalField() ) < small ){
        FatalErrorIn("createFieldRefs")
        << "Hydrogen species mass fraction cannot be zero anywhere in the domain" 
        << endl << exit(FatalError);
    }

    // Create lists with stoichiometric coefficients of species for reaction enthalpy
    for( label i=0; i<chemistry.nReaction(); i++ )
    {
        products_rrn.set
        (
            i,
            new List<specieCoeffs>(reactions[i].rhs())
        );
        reactants_rrn.set
        (
            i,
            new List<specieCoeffs>(reactions[i].lhs())
        );
    }

    // Compute enthalpy of reaction at standard conditions for each reaction
    for( label i=0; i<chemistry.nReaction(); i++ )
    {
        // Calculate enthalpy of reaction at standard conditions
        // This is Delta H_0R = H_{f0, products} - H_{f0, reactants}
        // Note the resulting units in J/kmol
        forAll(products_rrn[i], spk)
        {
            H0R[i] += products_rrn[i][spk].stoichCoeff*composition.Hf( products_rrn[i][spk].index )*composition.Wi(products_rrn[i][spk].index);
        }
        forAll(reactants_rrn[i], spk)
        {
            H0R[i] -= reactants_rrn[i][spk].stoichCoeff*composition.Hf( reactants_rrn[i][spk].index )*composition.Wi(reactants_rrn[i][spk].index);
        }
    }
}