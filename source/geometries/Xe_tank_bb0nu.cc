// ----------------------------------------------------------------------------
// nexus | Xe_tank_bb0nu.cc
//
// Tank filled with xenon for bb0nu decay simulation.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Xe_tank_bb0nu.h"

#include "CylinderPointSampler2020.h"
#include "MaterialsList.h"
#include "Visibilities.h"
#include "OpticalMaterialProperties.h"
#include "IonizationSD.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>
#include <G4VUserDetectorConstruction.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(Xe_tank_bb0nu, GeometryBase)

namespace nexus {

  Xe_tank_bb0nu::Xe_tank_bb0nu():
    GeometryBase(),
    liquid_(false),
    pressure_(STP_Pressure),
    radius_(1.*m),
    length_(3.*m)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/Xe_tank_bb0nu/",
      "Control commands of geometry Xenon tank for bb0nu simulations.");

    msg_->DeclareProperty("LXe", liquid_,
      "Build the tank with liquid xenon.");

    G4GenericMessenger::Command& pressure_cmd =
      msg_->DeclareProperty("pressure", pressure_,
      "Set pressure for gaseous xenon (if selected).");
    pressure_cmd.SetUnitCategory("Pressure");
    pressure_cmd.SetParameterName("pressure", false);
    pressure_cmd.SetRange("pressure>0.");

    G4GenericMessenger::Command& radius_cmd =
      msg_->DeclareProperty("radius", radius_, "Radius of the xenon tank.");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius", false);
    radius_cmd.SetRange("radius>0.");

    G4GenericMessenger::Command& length_cmd =
      msg_->DeclareProperty("length", length_, "Length of the xenon tank.");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length", false);
    length_cmd.SetRange("length>0.");

  }



  Xe_tank_bb0nu::~Xe_tank_bb0nu()
  {
    delete msg_;
  }



  void Xe_tank_bb0nu::Construct()
  {

    // INFO COUT________________________________________________________________

    std::cout<<"Tank size = "<<radius_<<" (radius) x "<<length_<<" (length)"<<std::endl;

    // LAB CREATION_____________________________________________________________

    // G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);
    G4double lab_z_ = length_ * 2;
    G4double lab_xy_ = radius_ * 2;
    G4Box* lab_solid = new G4Box("LAB", lab_xy_, lab_xy_, lab_z_);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);

    // Create a vertex generator for the tank
    tank_vertex_gen_ = new CylinderPointSampler2020(0., radius_, length_/2., 0, 2 * pi);


    // GENERAL PARAMETERS______________________________________________________

    G4double x = 0.;
    G4double y = 0.;
    G4double z = 0.;
    G4double rot_angle_;


    // Xe tank__________________________________________________________________

    G4String name = "Xe_tank";

    // Define solid volume as a cylinder
    G4Tubs* tank_solid =
     new G4Tubs(name, 0., radius_, length_/2., 0., 360.*deg);

    // Define the material (LXe or GXe) for the tank.
    // We use for this the NIST manager or the nexus materials list.
    G4Material* xenon = 0;
    if (liquid_)
      xenon = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");
    else
      xenon = materials::GXe(pressure_);


    // Define the logical volume of the tank using the material
    // and the solid volume defined above
    G4LogicalVolume* tank_logic =
    new G4LogicalVolume(tank_solid, xenon, name);
    GeometryBase::SetLogicalVolume(tank_logic);

    G4VisAttributes tank_col = nexus::Red();
    tank_logic->SetVisAttributes(tank_col);

    G4OpticalSurface* tank_opsur =
      new G4OpticalSurface("tank_opsur", unified, polished, dielectric_metal);
      // opsur->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());
      tank_opsur->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("tank_opsur", tank_logic, tank_opsur);

    // Set the logical volume of the tank as an ionization
    // sensitive detector, i.e. position, time and energy deposition
    // will be stored for each step of any charged particle crossing
    // the volume.
    IonizationSD* ionizsd = new IonizationSD("/Xe_tank");
    G4SDManager::GetSDMpointer()->AddNewDetector(ionizsd);
    tank_logic->SetSensitiveDetector(ionizsd);


    G4ThreeVector tank_pos = G4ThreeVector(x, y, z);

    G4RotationMatrix* tank_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    rot_angle_ = 0.;
    tank_rot_->rotateY(rot_angle_);

    new G4PVPlacement(G4Transform3D(*tank_rot_, tank_pos),
                      tank_logic, name, lab_logic,
                      false, 0, false);

    std::cout<<"Tank position = "<<tank_pos<<std::endl;


    // Al disk__________________________________________________________________

    G4String disk_name = "Al_disk";

    G4Material* disk_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    disk_mat->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    G4double disk_thickness = 1. * m;

    G4Tubs* disk_solid_vol =
      new G4Tubs(disk_name, 0., radius_, disk_thickness/2., 0., 360.*deg);

    G4LogicalVolume* disk_logic_vol =
      new G4LogicalVolume(disk_solid_vol, disk_mat, disk_name);

    G4OpticalSurface* Al_opsur =
      new G4OpticalSurface("Al_opsurf", unified, polished, dielectric_metal);
      // opsur->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());
      Al_opsur->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("Al_opsurf", disk_logic_vol, Al_opsur);

    // G4VisAttributes disk_col = nexus::LightBlue();
    G4VisAttributes disk_col = nexus::Blue();
    disk_logic_vol->SetVisAttributes(disk_col);
    // disk_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

    // G4ThreeVector disk_pos = G4ThreeVector(x - length_/2 - disk_thickness/2., y, z);
    G4ThreeVector disk_pos = G4ThreeVector(x, y, z);
    G4RotationMatrix* disk_rot_ = new G4RotationMatrix();
    rot_angle_ = pi/2.;
    // rot_angle_ = 0.;
    disk_rot_->rotateY(rot_angle_);
    new G4PVPlacement(G4Transform3D(*disk_rot_, disk_pos),
                      disk_logic_vol, disk_name, lab_logic,
                      false, 0, false);

    std::cout<<"AlDisk position = "<<disk_pos<<std::endl;



  }



  G4ThreeVector Xe_tank_bb0nu::GenerateVertex(const G4String& region) const
  {
    return tank_vertex_gen_->GenerateVertex(region);
  }


} // end namespace nexus
