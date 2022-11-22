// we'll try to work with more than one fiber

#include "Fiber_cylinder.h"
#include "GenericWLSFiber.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
// #include "CylinderPointSampler.h"
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

REGISTER_CLASS(Fiber_cylinder,GeometryBase)

// Fiber_cylinder::Fiber_cylinder():GeometryBase(), radius_(1.*mm), length_(1.*cm), cyl_vertex_gen_(0) {
Fiber_cylinder::Fiber_cylinder():GeometryBase(), radius_(1.*mm), length_(1.*cm),  n_fibers_(10){
    msg_=new G4GenericMessenger(this,"/Geometry/Fiber_cylinder/","Control commands of geometry OpticalFibre");

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

    msg_->DeclareProperty("n_fibers",n_fibers_,"Number of optical fibres");

    // cyl_vertex_gen_ = new CylinderPointSampler(0.5*radius_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);
}
Fiber_cylinder::~Fiber_cylinder() {
    // delete cyl_vertex_gen_;
    delete msg_;
}
void Fiber_cylinder::Construct(){
    // G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);
    G4Box* lab_solid = new G4Box("LAB", 20 * mm,20 * mm,10.1*cm);

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

    // fibers loop
    // int N = 70;
    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;
    // G4double dif_phi = .6*2*pi/n_fibers_; // angular separation between fibers
    G4double dif_phi = .5*2*pi/n_fibers_; // angular separation between fibers
    // G4double dif_phi = 2*pi/n_fibers_; // angular separation between fibers
    G4double rad = n_fibers_*radius_/pi; // radii of the big cylinder of fibers
    G4double ang_pos;
    G4double x;
    G4double y;
    G4double z = 0.;

    // for (int i=0; i<=2*n_fibers_; i++){
    for (int i=0; i<=n_fibers_/0.5; i++){
      // fiber3
      fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
      fiber_->SetCoreOpticalProperties(opticalprops::Y11());
      fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
      fiber_->Construct();
      fiber_logic = fiber_->GetLogicalVolume();

      ang_pos = dif_phi*i;
      x = rad * cos(ang_pos);
      y = rad * sin(ang_pos);
      new G4PVPlacement(0,G4ThreeVector(x, y, z),fiber_logic,
                              fiber_logic->GetName(),lab_logic,true,0,true);
    }


}


G4ThreeVector Fiber_cylinder::GenerateVertex(const G4String& region) const {
    // return cyl_vertex_gen_->GenerateVertex(region);
    // G4ThreeVector vertex(1.,1.,1.);
    G4ThreeVector vertex(0., 0., 0.);

    // WORLD
    if (region == "WHOLE_VOL") {
      return vertex;
    }
    // else if (region == "AD_HOC") {
    //   return specific_vertex_;
    // }
    else {
      G4Exception("[Fiber_cylinder]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }
    return vertex;

}
