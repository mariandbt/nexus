// ----------------------------------------------------------------------------
// nexus | Fib_box_struct.cc
//
// Box structure to PMTs calibration using WSL fibers
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Fib_box_struct.h"
#include "GenericWLSFiber.h"
#include "SiPM11_eff.h"

#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4Colour.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4NistManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(Fib_box_struct,GeometryBase)

Fib_box_struct::Fib_box_struct():
    GeometryBase(),
    radius_(1.*mm),
    length_(1.*cm),
    box_xy_(36.*mm),
    box_z_(14.*cm)
  {
    std::cout<<"HERE!"<<std::endl;
    msg_=new G4GenericMessenger(this,"/Geometry/Fib_box_struct/",
        "Control commands of geometry Fiber box structure.");

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

    // G4GenericMessenger::Command& box_xy_cmd =
    //         msg_->DeclareProperty("box_xy",box_xy_,"Side of the fiber box structure");
    // radius_cmd.SetUnitCategory("Length");
    // radius_cmd.SetParameterName("box_xy",false);
    // radius_cmd.SetRange("box_xy>0.");
    //
    // G4GenericMessenger::Command& box_z_cmd =
    //         msg_->DeclareProperty("box_z",box_z_,"Length of the fiber box structure");
    // length_cmd.SetUnitCategory("Length");
    // length_cmd.SetParameterName("box_z",false);
    // length_cmd.SetRange("box_z>0.");


    sipm_ = new SiPM11_eff();

}
Fib_box_struct::~Fib_box_struct() {
    delete msg_;
}
void Fib_box_struct::Construct(){
    // LAB CREATION___________________________________________________

    // G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);
    G4double lab_z_ = box_z_ * 2;
    G4double lab_xy_ = length_ * 2;
    G4Box* lab_solid = new G4Box("LAB", lab_xy_, lab_xy_, lab_z_);

    // cyl_vertex_gen_ = new CylinderPointSampler2020(0., radius_cyl_, length_/2, 0, 2 * pi);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);

    G4double x = 0.;
    G4double y = 0.;
    G4double z = 0.;


    // fibers______________________________________________________

    // G4Material* ps = materials::PS();
    G4Material* ps = materials::Y11();
    G4Material* tpb = materials::TPB();

    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;
    // the factor 2 it's because "radius_" actually stands for the fiber DIAMETER
    G4double n_fibers = 33;
    // std::cout<<"n_fibers = "<<n_fibers<<std::endl;

    fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
    fiber_->SetCoreOpticalProperties(opticalprops::Y11());
    fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
    fiber_->Construct();
    fiber_logic = fiber_->GetLogicalVolume();


    // SiPM______________________________________________________

    G4LogicalVolume* sipm_logic;

    sipm_->Construct();
    sipm_logic = sipm_->GetLogicalVolume();


    // Al DISK____________________________________________________

    G4String disk_name = "Al DISK";

    G4Material* disk_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    disk_mat->SetMaterialPropertiesTable(opticalprops::PolishedAl());


    G4double disk_thickness = .1 * mm;

    G4Tubs* disk_solid_vol =
      new G4Tubs(disk_name, 0., radius_/2., disk_thickness/2., 0., 360.*deg);

    G4LogicalVolume* disk_logic_vol =
      new G4LogicalVolume(disk_solid_vol, disk_mat, disk_name);

    G4OpticalSurface* disk_opsur =
      new G4OpticalSurface("Al_OPSURF", unified, polished, dielectric_metal);
      // disk_opsur->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());
      disk_opsur->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("Al_OPSURF", disk_logic_vol, disk_opsur);

    // G4VisAttributes disk_col = nexus::LightBlue();
    G4VisAttributes disk_col = nexus::Blue();
    disk_logic_vol->SetVisAttributes(disk_col);
    // disk_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());


    // Black box_____________________________________________________

    G4String box_name = "Black box";
    G4String box_side_name = "Black box side";

    G4Box* full_box_solid_vol =
      new G4Box(box_name, box_xy_/2., box_xy_/2., box_z_/2.);

    G4double side_thickness = .1 * mm;
    G4Box* box_side_solid_vol =
    // new G4Box(box_side_name, box_xy_/2., box_xy_/2., side_thickness/2.);
    new G4Box(box_side_name, box_xy_/2., box_xy_/4., box_z_/4.);

    G4ThreeVector side_pos = G4ThreeVector(0., 0., box_z_/2. + side_thickness/2.);
    G4RotationMatrix* side_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    G4double rot_angle_ = 0.;
    side_rot_->rotateY(rot_angle_);
    G4SubtractionSolid* box_solid_vol =
      new G4SubtractionSolid("Box-Side", full_box_solid_vol, box_side_solid_vol,
                              side_rot_, side_pos);

    G4Material* box_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");
    // box_mat->SetMaterialPropertiesTable(opticalprops::PTFE());

    G4LogicalVolume* box_logic_vol =
      new G4LogicalVolume(box_solid_vol, box_mat, box_name);
    G4VisAttributes box_col = nexus::Red();
    // box_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());
    box_logic_vol->SetVisAttributes(box_col);

    G4ThreeVector box_pos = G4ThreeVector(0., 0., box_z_/2.);

    G4RotationMatrix* box_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    rot_angle_ = 0.;
    box_rot_->rotateY(rot_angle_);

    new G4PVPlacement(G4Transform3D(*box_rot_, box_pos),
                      box_logic_vol, box_name, lab_logic,
                      false, 0, false);

    // Teflon panel_____________________________________________________

    G4String panel_name = "Teflon panel";

    G4double panel_thickness = .1 * mm;

    G4Material* panel_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    panel_mat->SetMaterialPropertiesTable(opticalprops::PTFE());

    G4Box* panel_solid_vol =
      new G4Box(panel_name, box_xy_/2., box_xy_/2., panel_thickness/2.);

    G4LogicalVolume* panel_logic_vol =
      new G4LogicalVolume(panel_solid_vol, panel_mat, panel_name);
    G4VisAttributes panel_col = nexus::White();
    panel_logic_vol->SetVisAttributes(panel_col);

    G4OpticalSurface* panel_opsur =
      new G4OpticalSurface("Al_OPSURF", unified, polished, dielectric_metal);
      // disk_opsur->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());
      panel_opsur->SetMaterialPropertiesTable(opticalprops::PTFE());
    new G4LogicalSkinSurface("Al_OPSURF", panel_logic_vol, panel_opsur);

    G4double panel_z_pos = box_z_ + radius_ + panel_thickness/2.;
    G4ThreeVector panel_pos = G4ThreeVector(0., 0., panel_z_pos);

    G4RotationMatrix* panel_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    rot_angle_ = 0.;
    panel_rot_->rotateY(rot_angle_);
    // new G4PVPlacement(G4Transform3D(*panel_rot_, panel_pos),
    //                   panel_logic_vol, panel_name, lab_logic,
    //                   false, 0, false);

    // loop_____________________________________________________________

    for (int i=0; i < n_fibers; i++){

      // // fiber
      // fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
      // fiber_->SetCoreOpticalProperties(opticalprops::Y11());
      // fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
      // fiber_->Construct();
      // fiber_logic = fiber_->GetLogicalVolume();

      y = i*radius_ - box_xy_/2;
      z = box_z_ + radius_/2;

      G4RotationMatrix* fib_rot_ = new G4RotationMatrix();
      rot_angle_ = pi/2.;
      // rot_angle_ = 0.;
      fib_rot_->rotateY(rot_angle_);

      // new G4PVPlacement(fib_rot_,G4ThreeVector(x, y, z),fiber_logic,
      //                         fiber_logic->GetName(),lab_logic,true,0,true);

      // // SiPM
      // sipm_->Construct();
      // sipm_logic = sipm_->GetLogicalVolume();

      // to avoid overlap among SiPMs intercalate them in X
      G4double sipm_x_pos = x + length_/2. + .45 * mm;
      // std::cout<<"sipm_x_pos = "<<sipm_x_pos<<std::endl;
      G4ThreeVector sipm_pos = G4ThreeVector(sipm_x_pos, y, z);

      G4RotationMatrix* sipm_rot_ = new G4RotationMatrix();
      rot_angle_ = pi;
      // rot_angle_ = 0.;
      sipm_rot_->rotateY(rot_angle_);
      // new G4PVPlacement(G4Transform3D(*sipm_rot_, sipm_pos), sipm_logic,
      //                   sipm_logic->GetName(),lab_logic,true,0,true);


      // Al disk

      G4ThreeVector disk_pos = G4ThreeVector(x - length_/2 - disk_thickness/2., y, z);

      // new G4PVPlacement(0, disk_pos,
      //                   disk_logic_vol, disk_name, lab_logic,
      //                   false, 0, false);

    }

    // std::cout<<"R = "<<radius_cyl_<<std::endl;
    // std::cout<<"L = "<<length_<<std::endl;

}
G4ThreeVector Fib_box_struct::GenerateVertex(const G4String& region) const {
    // return cyl_vertex_gen_->GenerateVertex(region);

    // // G4ThreeVector vertex(1.,1.,1.);
    // G4ThreeVector vertex(box_xy_/2, box_xy_/2, 0.);
    G4ThreeVector vertex(0., 0., 0.);

    // WORLD
    if (region == "CENTER") {
      return vertex;
    }

    return vertex;

}
