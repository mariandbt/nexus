#include "Fib_SiPM_cyl.h"
#include "GenericWLSFiber.h"
#include "SiPM11_eff.h"

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

REGISTER_CLASS(Fib_SiPM_cyl,GeometryBase)

Fib_SiPM_cyl::Fib_SiPM_cyl():
    GeometryBase(),
    radius_(1.*mm),
    length_(1.*cm),
    radius_cyl_(1. *cm),
    // sipm_visibility_(1),
    cyl_vertex_gen_(0)
  {
    msg_=new G4GenericMessenger(this,"/Geometry/Fib_SiPM_cyl/",
        "Control commands of geometry Fiber cylinder with SiPMs.");

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

    G4GenericMessenger::Command& radius_cyl_cmd =
            msg_->DeclareProperty("radius_cyl",radius_cyl_,"Radius of the cylinder");
    radius_cyl_cmd.SetUnitCategory("Length");
    radius_cyl_cmd.SetParameterName("radius_cyl",false);
    radius_cyl_cmd.SetRange("radius_cyl>0.");

    cyl_vertex_gen_ = new CylinderPointSampler(0.5*radius_cyl_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);

    sipm_ = new SiPM11_eff();

}
Fib_SiPM_cyl::~Fib_SiPM_cyl() {
    // delete sipm_;
    delete cyl_vertex_gen_;
    delete msg_;
}
void Fib_SiPM_cyl::Construct(){
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

    // fibers
    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;
    G4double n_fibers = radius_cyl_*pi/radius_;
    G4double dif_theta = .5*2*pi/n_fibers; // angular separation between fibers
    G4double theta;
    G4double x;
    G4double y;
    G4double z = 0.;

    // SiPM
    G4LogicalVolume* sipm_logic;

    // loop
    for (int i=0; i<=2*n_fibers; i++){

      // fiber
      fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
      fiber_->SetCoreOpticalProperties(opticalprops::Y11());
      fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
      fiber_->Construct();
      fiber_logic = fiber_->GetLogicalVolume();

      theta = dif_theta*i;
      x = radius_cyl_ * cos(theta);
      y = radius_cyl_ * sin(theta);
      new G4PVPlacement(0,G4ThreeVector(x, y, z),fiber_logic,
                              fiber_logic->GetName(),lab_logic,true,0,true);

      // SiPM
      sipm_->Construct();
      sipm_logic = sipm_->GetLogicalVolume();

      G4ThreeVector sipm_pos = G4ThreeVector(x, y, z + length_/2);

      G4RotationMatrix* sipm_rot_ = new G4RotationMatrix();
      G4double rot_angle_ = pi;
      //G4double rot_angle_ = 0.;
      sipm_rot_->rotateY(rot_angle_);
      new G4PVPlacement(G4Transform3D(*sipm_rot_, sipm_pos), sipm_logic,
                        sipm_logic->GetName(),lab_logic,true,0,true);

    }


}


G4ThreeVector Fib_SiPM_cyl::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);

    // // G4ThreeVector vertex(1.,1.,1.);
    // G4ThreeVector vertex(0., 0., 0.);
    // return vertex;

}
