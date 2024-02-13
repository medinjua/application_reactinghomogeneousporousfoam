#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureReference pressureReference(p, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// ****************************************************************************
// Additions for porous catalyst simulations
// ****************************************************************************

// chemistryModel required
autoPtr<basicChemistryModel> pChemistry(basicChemistryModel::New(thermo));

const word inertSpecie(thermo.properties().lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo.properties())
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

// Gas constant units in J/(kmol * K)
const dimensionedScalar RuSI  = dimensionedScalar(
                                    dimensionSet(1,2,-2,-1,-1,0,0), 
                                    8314.463
                                );

Info<< "\nAdding porosity" << endl;
volScalarField epspore
(
    IOobject
    (
        "epspore",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "\nAdding effective thermal conductivity" << endl;
volTensorField kappaeff
(
    IOobject
    (
        "kappaeff",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "\nAdding solid density and specific heat capacity" << endl;
volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField cpvs
(
    IOobject
    (
        "cpvs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "\nAdding catalyst density (dimensionless, see rhocat file)" << endl;
volScalarField rhocat
(
    IOobject
    (
        "rhocat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Inhibition parameter field
volScalarField DENfact
(
    IOobject
    (
        "DENfact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DENfact",dimensionSet(0,0,0,0,0,0, 0),0)
); 
// NOTE: Concentration DEN factor is not dimensionless
// It is [kmol/m3]^n, where n is the exponent indicated by expH2DENConc
// However, for ease of use of the function pow(a,b) = a^b, 
// DENConc is assumed dimless
volScalarField DENConc
(
    IOobject
    (
        "DENConc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DENConc",dimensionSet(0,0,0,0,0,0, 0),0)
); 

// Heat release rate [J/m3/s]
volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Field for species i reaction rate [kg/m3/s]
volScalarField omegarr_i
(
    IOobject
    (
        "omegarr_i",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("omegarr_i",dimensionSet(1,-3,-1,0,0,0, 0),0)
);