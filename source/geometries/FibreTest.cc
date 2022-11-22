#include "FibreTest.h"
#include "GenericWLSFiber.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
#include "CylinderPointSampler.h"
#include "FactoryBase.h"

#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4NistManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(FibreTest,GeometryBase)

FibreTest::FibreTest():GeometryBase(), radius_(1.*mm), length_(1.*cm), cyl_vertex_gen_(0) {
    msg_=new G4GenericMessenger(this,"/Geometry/FibreTest/","Control commands of geometry OpticalFibre");

    G4GenericMessenger::Command& radius_cmd = 
            msg_->DeclareProperty("radius",radius_,"Radius of the cylindrical optical fibre");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius",false);
    radius_cmd.SetRange("radius>0.");

    G4GenericMessenger::Command& length_cmd = 
            msg_->DeclareProperty("length",length_,"Length of the cylindrical optical fibre");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length",false);
    length_cmd.SetRange("length>0.");


    cyl_vertex_gen_ = new CylinderPointSampler(0.5*radius_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);
}
FibreTest::~FibreTest() {
    delete cyl_vertex_gen_;
    delete msg_;
}
void FibreTest::Construct(){
    G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);
    G4Material* ps = materials::PS();
    G4Material* tpb = materials::TPB();
    GenericWLSFiber* fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
    fiber_->SetCoreOpticalProperties(opticalprops::Y11());
    fiber_->Construct();
    G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
    new G4PVPlacement(0,G4ThreeVector(0,0,0),fiber_logic,
                            fiber_logic->GetName(),lab_logic,true,0,true);
}

G4ThreeVector FibreTest::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);
}