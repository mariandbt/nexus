// ----------------------------------------------------------------------------
// nexus | FibBarrMeth.cc
//
// Cylinder containing optical fibers with a methacrylate window
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------


#include "FibBarrMeth.h"

#include "Visibilities.h"
#include "FactoryBase.h"
#include "OpticalMaterialProperties.h"

#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Colour.hh>
#include <G4Material.hh>


#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>

using namespace nexus;

REGISTER_CLASS(FibBarrMeth,GeometryBase)

namespace nexus {

FibBarrMeth::FibBarrMeth():
GeometryBase(),
fiber_type_ ("Y11"),
coated_ (true),
diameter_(1.*mm),
length_(1.*cm),
radius_cyl_(1. *cm)//,
// sensor_visibility_ (true)
{
    msg_=new G4GenericMessenger(this,"/Geometry/FibBarrMeth/",
                                "Control commands of geometry FibBarrMeth");

    msg_->DeclareProperty("fiber_type", fiber_type_,
            "Fiber type");

    msg_->DeclareProperty("coated", coated_,
          "Option to coat or not the fibers");

    G4GenericMessenger::Command& diameter_cmd =
            msg_->DeclareProperty("diameter",diameter_,"Diameter of the cylindrical optical fibre");
    diameter_cmd.SetUnitCategory("Length");
    diameter_cmd.SetParameterName("diameter",false);
    diameter_cmd.SetRange("diameter>0.");

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
    //
    // msg_->DeclareProperty("sensor_visibility", sensor_visibility_,
    //                       "Sensors visibility");

}
FibBarrMeth::~FibBarrMeth() {
    delete msg_;
}
void FibBarrMeth::Construct(){

    // LAB CREATION___________________________________________________

    cyl_vertex_gen_ = new CylinderPointSampler2020(0., radius_cyl_, length_/2, 0, 2*pi);

    G4double teflon_thickness_ = 5. * mm;

    G4double lab_z_ = (radius_cyl_ + diameter_/2. + teflon_thickness_)*2 ;
    G4double lab_xy_ = (length_ + 2*teflon_thickness_ + 2*mm) * 2;
    G4Box* lab_solid = new G4Box("LAB", lab_xy_, lab_xy_, lab_z_);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    // this->SetLogicalVolume(lab_logic);
    GeometryBase::SetLogicalVolume(lab_logic);


    G4double x;
    G4double y;
    G4double z = 0.;
    G4double rot_angle;

    G4Material *this_coating = nullptr;
    G4MaterialPropertiesTable *this_coating_optical = nullptr;

    G4Material *this_fiber = nullptr;
    G4MaterialPropertiesTable *this_fiber_optical = nullptr;

    if (fiber_type_ == "Y11") {
      this_fiber = materials::Y11();
      this_fiber_optical = opticalprops::Y11();

      this_coating = materials::TPB();
      this_coating_optical = opticalprops::TPB();

    } else if (fiber_type_ == "B2") {

      this_fiber = materials::B2();
      this_fiber_optical = opticalprops::B2();

      this_coating = materials::TPH();
      this_coating_optical = opticalprops::TPH();

    } else {
      G4Exception("[FiberBarrel]", "Construct()",
                  FatalException, "Invalid fiber type, must be Y11 or B2");
    }



    // Al disk _____________________________________________________________________

   G4String disk_name = "Al_DISK";

   G4Material* disk_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  //  disk_mat->SetMaterialPropertiesTable(opticalprops::PolishedAl());

   G4double disk_thickness = .1 * mm;

   G4Tubs* disk_solid_vol =
     new G4Tubs(disk_name, 0., diameter_/2., disk_thickness/2., 0., 360.*deg);

   G4LogicalVolume* disk_logic_vol =
     new G4LogicalVolume(disk_solid_vol, disk_mat, disk_name);

   G4OpticalSurface* disk_opsur =
     new G4OpticalSurface("Al_OPSURF", unified, polished, dielectric_metal);
     disk_opsur->SetMaterialPropertiesTable(opticalprops::PolishedAl());

   new G4LogicalSkinSurface("Al_OPSURF", disk_logic_vol, disk_opsur);

   // G4VisAttributes disk_col = nexus::LightBlue();
   G4VisAttributes disk_col = nexus::Blue();
   disk_logic_vol->SetVisAttributes(disk_col);
   // disk_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());


   // FAKE DETECTOR /////////////////////////////////////////////

   G4Material* sens_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
   G4double sensor_thickness = 1. * mm;

   G4double sens_z = sensor_thickness;
  //  G4double sens_z = 0.1 * mm;

   G4Tubs* sens_solid_vol =
     new G4Tubs("disk", 0, diameter_/2, sens_z, 0, 2 * M_PI);

   G4LogicalVolume* sens_logic_vol =
     new G4LogicalVolume(sens_solid_vol, sens_mat, "DISK");

   G4OpticalSurface* opsur =
     new G4OpticalSurface("PERFECT_OPSURF", unified, polished, dielectric_metal);
   opsur->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());

   new G4LogicalSkinSurface("PERFECT_OPSURF", sens_logic_vol, opsur);


  //  // Sensor _____________________________________________________________________
  //
  //  G4double sensor_thickness = 1. * mm;
  //  G4String sensor_name = "F_SENSOR";
  //
  // /// Build the sensor
  // photo_sensor_  = new GenericCircularPhotosensor(sensor_name, diameter_, sensor_thickness);
  //
  //
  // /// Constructing the sensors
  // // Optical Properties of the sensor
  // G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
  //
  //   // perfect detector
  //   // G4int entries = 2;
  //   G4double energy[]       = {0.2 * eV, 11.5 * eV};
  //   G4double efficiency[]   = {1.      , 1.       };
  //   G4double reflectivity[] = {0.      , 0.       };
  //
  //   photosensor_mpt->AddProperty("EFFICIENCY", energy, efficiency, 2);
  //   photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 2);
  //
  //
  //   G4Material* fiber_mat = this_fiber;
  //   fiber_mat->SetMaterialPropertiesTable(this_fiber_optical);
  //
  //   G4MaterialPropertyVector* fibers_rindex =
  //   fiber_mat->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  //
  //   photo_sensor_ ->SetOpticalProperties(photosensor_mpt);
  //
  //   // Adding to sensors encasing, the Refractive Index of fibers to avoid reflections
  //   photo_sensor_ ->SetWindowRefractiveIndex(fibers_rindex);
  //
  //   // Setting the time binning
  //   // photo_sensor_ ->SetTimeBinning(100. * ns); // Size of fiber sensors time binning
  //
  //   // Set mother depth & naming order
  //   photo_sensor_ ->SetSensorDepth(1);
  //   photo_sensor_ ->SetMotherDepth(2);
  //   photo_sensor_ ->SetNamingOrder(1);
  //
  //   // Set visibilities
  //   photo_sensor_ ->SetVisibility(sensor_visibility_);
  //
  //   // Construct
  //   photo_sensor_ ->Construct();
  //
  //   G4LogicalVolume* photo_sensor_logic  = photo_sensor_ ->GetLogicalVolume();


    // Teflon caps______________________________________________________________________________

    G4String caps_name = "TEFLON_CAPS";

    G4double caps_thickness = teflon_thickness_;
    // G4double caps_thickness = sensor_thickness;

    G4Material* teflon_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    teflon_mat->SetMaterialPropertiesTable(opticalprops::PTFE());

    G4OpticalSurface* teflon_opsur =
    new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
    teflon_opsur->SetMaterialPropertiesTable(opticalprops::PTFE());


    G4Tubs* caps_solid_vol =
    new G4Tubs(caps_name, 0., radius_cyl_ - diameter_/2., caps_thickness/2., 0., 360.*deg);

    G4LogicalVolume* caps_logic_vol =
      new G4LogicalVolume(caps_solid_vol, teflon_mat, caps_name);

    new G4LogicalSkinSurface("TEFLON_CAPS_OPSURF", caps_logic_vol, teflon_opsur);

    G4VisAttributes teflon_col = nexus::White();
    caps_logic_vol->SetVisAttributes(teflon_col);
    // caps_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4RotationMatrix* caps_rot = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    caps_rot->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*caps_rot, G4ThreeVector(0., 0., (length_ + caps_thickness)/2.)),
                      caps_logic_vol, caps_name + "_1", lab_logic,
                      true, 0, false);
    new G4PVPlacement(G4Transform3D(*caps_rot, G4ThreeVector(0., 0., -(length_ + caps_thickness)/2.)),
                      caps_logic_vol, caps_name + "_2", lab_logic,
                      true, 1, false);


    // Outer teflon cilynder______________________________________________________________________________

    G4String teflon_cyl_name = "TEFLON_CYLINDER";

    G4Tubs* teflon_cyl_solid_vol =
    new G4Tubs(teflon_cyl_name, radius_cyl_ + diameter_/2., radius_cyl_ + diameter_/2. + teflon_thickness_,
      (length_ + sensor_thickness + disk_thickness)/2., 0., 360.*deg);

    G4LogicalVolume* teflon_cyl_logic_vol =
      new G4LogicalVolume(teflon_cyl_solid_vol, teflon_mat, teflon_cyl_name);

    new G4LogicalSkinSurface("TEFLON_CYLINDER_OPSURF", teflon_cyl_logic_vol, teflon_opsur);

    teflon_cyl_logic_vol->SetVisAttributes(teflon_col);
    // teflon_cyl_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4ThreeVector teflon_cyl_pos = G4ThreeVector(0., 0., 0.);
    G4RotationMatrix* teflon_cyl_rot = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    teflon_cyl_rot->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*teflon_cyl_rot, teflon_cyl_pos),
                      teflon_cyl_logic_vol, teflon_cyl_name, lab_logic,
                      false, 0, false);



    // Fiber _____________________________________________________________________

    G4double n_fibers = floor(2*radius_cyl_*pi/diameter_);
    G4double dif_theta = 2*pi/n_fibers; // angular separation between fibers
    G4double theta;
    G4cout<<"n_fibers = "<<n_fibers<<G4endl;


    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;

    fiber_ = new GenericWLSFiber(fiber_type_, true, diameter_, length_,
                                 true, coated_, this_coating, this_fiber, true);

    fiber_->SetCoreOpticalProperties(this_fiber_optical);
    if (coated_) {
      fiber_->SetCoatingOpticalProperties(this_coating_optical);
    }

    fiber_->Construct();
    fiber_logic = fiber_->GetLogicalVolume();


    // Loop to place elements

    for (G4int i=0; i < n_fibers; i++){

      theta = dif_theta*i;
      x = radius_cyl_ * cos(theta);
      y = radius_cyl_ * sin(theta);
      std::string label = std::to_string(i);

      // fibers
      new G4PVPlacement(0,G4ThreeVector(x, y, z),fiber_logic,
                        fiber_logic->GetName() + "_" + label,lab_logic,
                        true, i, false);

      // aluminization
      new G4PVPlacement(0,G4ThreeVector(x, y, z - length_/2. - disk_thickness/2.),
                        disk_logic_vol, disk_name + "_" + label,lab_logic,
                        true, n_fibers + i, false);

      // sensor
      new G4PVPlacement(0, G4ThreeVector(x, y, length_/2 + sens_z),
                        sens_logic_vol, "DISKL-" + label, lab_logic,
                        false, n_fibers+i, false);

      // G4double sensor_z_pos =  z + length_/2. + sensor_thickness/2.;
      // G4ThreeVector sensor_pos = G4ThreeVector(x, y, sensor_z_pos);
      //
      // G4RotationMatrix* sensor_rot = new G4RotationMatrix();
      // // rot_angle = -pi/2.;
      // // rot_angle = 0.;
      // rot_angle = pi;
      // sensor_rot->rotateY(rot_angle);
      // new G4PVPlacement(G4Transform3D(*sensor_rot, sensor_pos), photo_sensor_logic,
      //                   photo_sensor_logic->GetName() + "_" + label,lab_logic,
      //                   true, 2*n_fibers + i, false);



    }


}


G4ThreeVector FibBarrMeth::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);
  }
}
