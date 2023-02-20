// ----------------------------------------------------------------------------
// nexus | Fib_box_struct.cc
//
// Box structure to PMTs calibration using WSL fibers
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Fib_box_struct.h"
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"

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
    box_xy_(40.*mm),
    box_z_(14.*cm),
    side_thickness_ (2. * mm),
    sensor_type_ ("PERFECT"),
    sensor_size_ (40. * mm),
    sensor_thickness_ (2. * mm),
    sensor_visibility_ (true)
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

    G4GenericMessenger::Command& sensor_size_cmd =
            msg_->DeclareProperty("sensor_size",sensor_size_,"Side length of squared fiber sensors");
    sensor_size_cmd.SetUnitCategory("Length");
    sensor_size_cmd.SetParameterName("sensor_size",false);
    sensor_size_cmd.SetRange("sensor_size>0.");

    msg_->DeclareProperty("sensor_visibility", sensor_visibility_,
                          "Sensors visibility");

    msg_->DeclareProperty("sensor_type", sensor_type_,
                          "Sensors type");

}
Fib_box_struct::~Fib_box_struct() {
    delete msg_;
}
void Fib_box_struct::Construct(){

    // INTRO COUT___________________________________________________

    std::cout<<"Selected sensor type = "<<sensor_type_<<std::endl;
    std::cout<<"Sensor size = "<<sensor_size_<<std::endl;


    // LAB CREATION___________________________________________________

    // G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);
    G4double lab_z_ = box_z_ * 2;
    G4double lab_xy_ = length_ * 2;
    G4Box* lab_solid = new G4Box("LAB", lab_xy_, lab_xy_, lab_z_);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);


    // GENERAL PARAMETERS______________________________________________________

    G4double x = length_/2 - box_xy_/2 - 1*mm; // x-position of the fibers
    // G4double x = 0.; // x position of the fibers
    G4double y = 0.;
    G4double z = 0.;
    G4double rot_angle_;


    // FIBERS______________________________________________________

    G4Material* ps = materials::Y11();
    G4Material* tpb = materials::TPB();

    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;
    // IMPORTANT: "radius_" actually stands for the fiber DIAMETER
    G4double n_fibers = 33;
    // std::cout<<"n_fibers = "<<n_fibers<<std::endl;

    fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
    fiber_->SetCoreOpticalProperties(opticalprops::Y11());
    // fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
    fiber_->Construct();
    fiber_logic = fiber_->GetLogicalVolume();


    // OPTICAL GEL____________________________________________________

    G4String opt_gel_name = "Optical gel";

    G4Material* optical_coupler = materials::OpticalSilicone();
    optical_coupler->SetMaterialPropertiesTable(opticalprops::OptCoupler());

    G4double opt_gel_thickness = .1 * mm;

    G4Tubs* opt_gel_solid_vol =
      new G4Tubs(opt_gel_name, 0., radius_/2., opt_gel_thickness/2., 0., 360.*deg);

    G4LogicalVolume* opt_gel_logic_vol =
      new G4LogicalVolume(opt_gel_solid_vol, optical_coupler, opt_gel_name);

    G4VisAttributes opt_gel_col = nexus::Lilla();
    opt_gel_logic_vol->SetVisAttributes(opt_gel_col);
    // opt_gel_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());


    // SENSOR______________________________________________________

    /// Build the sensor
    photo_sensor_  = new GenericPhotosensor("F_SENSOR", sensor_size_,
                                            sensor_size_, sensor_thickness_);
    /// Constructing the sensors
    // Optical Properties of the sensor
    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();

    const G4int entries = 9;
    G4double energy[entries];
    G4double reflectivity[entries];
    G4double efficiency[entries];

    if (sensor_type_ == "PERFECT") {
      // perfect detector
      G4double energy_values[entries]       = {0.2 * eV, 3.5 * eV, 3.6 * eV,
                                               4.2 * eV, 5.5 * eV, 7.6 * eV,
                                               8.2 * eV, 9.5 * eV, 11.5 * eV};
      G4double reflectivity_values[entries] = {0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      };
      G4double efficiency_values[entries]   = {1.      , 1.      , 1.      ,
                                               1.      , 1.      , 1.      ,
                                               1.      , 1.      , 1.      };

      for (G4int n=0; n<entries; ++n) {
        energy[n]       = energy_values[n];
        reflectivity[n] = reflectivity_values[n];
        efficiency[n]   = efficiency_values[n];
      }
    }

    else if (sensor_type_ == "SiPM") {
      // SiPM
      G4double energy_values[entries]       = {
        h_Planck * c_light / (815.514 * nm), h_Planck * c_light / (725.221 * nm),
        h_Planck * c_light / (626.346 * nm), h_Planck * c_light / (545.181 * nm),
        h_Planck * c_light / (457.196 * nm), h_Planck * c_light / (392.546 * nm),
        h_Planck * c_light / (369.434 * nm), h_Planck * c_light / (343.040 * nm),
        h_Planck * c_light / (328.480 * nm)
      };
      G4double reflectivity_values[entries] = {0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      };
      G4double efficiency_values[entries]   = {.09049, .16291, .28700, .41550,
                                               .49751, .43264, .34189, .20612,
                                               .06887
                                              };

      for (G4int n=0; n<entries; ++n) {
        energy[n]       = energy_values[n];
        reflectivity[n] = reflectivity_values[n];
        efficiency[n]   = efficiency_values[n];
      }
    }

    else if (sensor_type_ == "PMT") {
      // PMT
      G4double energy_values[entries]       = {
        h_Planck * c_light / (903.715 * nm), h_Planck * c_light / (895.975 * nm),
        h_Planck * c_light / (866.563 * nm), h_Planck * c_light / (826.316 * nm),
        h_Planck * c_light / (628.173 * nm), h_Planck * c_light / (490.402 * nm),
        h_Planck * c_light / (389.783 * nm), h_Planck * c_light / (330.96 * nm),
        h_Planck * c_light / (296.904 * nm)
      };
      G4double reflectivity_values[entries] = {0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      ,
                                               0.      , 0.      , 0.      };
      G4double efficiency_values[entries]   = {0.00041, 0.00107, 0.01248, 0.06181,
                                               0.12887, 0.19246, 0.09477, 0.06040,
                                               0.00826};

      for (G4int n=0; n<entries; ++n) {
        energy[n]       = energy_values[n];
        reflectivity[n] = reflectivity_values[n];
        efficiency[n]   = efficiency_values[n];
      }
    }

    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, entries);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   entries);
    photo_sensor_ ->SetOpticalProperties(photosensor_mpt);

    // Adding to sensors encasing, the Refractive Index of fibers to avoid reflections
    G4MaterialPropertyVector* fibers_rindex =
      ps->GetMaterialPropertiesTable()->GetProperty("RINDEX");
    photo_sensor_ ->SetWindowRefractiveIndex(fibers_rindex);

    // Setting the time binning
    // photo_sensor_ ->SetTimeBinning(100. * ns); // Size of fiber sensors time binning

    // Set mother depth & naming order
    photo_sensor_ ->SetSensorDepth(1);
    photo_sensor_ ->SetMotherDepth(2);
    photo_sensor_ ->SetNamingOrder(1);

    // Set visibilities
    photo_sensor_ ->SetVisibility(sensor_visibility_);

    // Construct
    photo_sensor_ ->Construct();

    G4LogicalVolume* photo_sensor_logic  = photo_sensor_ ->GetLogicalVolume();

    // Sensor placement
    G4double sensor_x_pos = x + length_/2. + opt_gel_thickness/2. + sensor_thickness_/2.;

    if (sensor_type_ == "SiPM") {

      // Metacrilate window

      G4String window_name = "Metacrilate window";

      G4double window_thickness = 1.5 * mm;

      G4Material* window_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
      window_mat->SetMaterialPropertiesTable(opticalprops::GlassEpoxy());

      G4Box* window_solid_vol =
        new G4Box(window_name, box_xy_/2., box_xy_/2., window_thickness/2.);

      G4LogicalVolume* window_logic_vol =
        new G4LogicalVolume(window_solid_vol, window_mat, window_name);
      G4VisAttributes window_col = nexus::LightBlue();
      window_logic_vol->SetVisAttributes(window_col);

      G4double window_x_pos = x + length_/2. + opt_gel_thickness/2. + window_thickness/2.;
      G4ThreeVector window_pos = G4ThreeVector(window_x_pos, 0., box_z_);

      G4RotationMatrix* window_rot_ = new G4RotationMatrix();
      rot_angle_ = pi/2.;
      // rot_angle_ = 0.;
      window_rot_->rotateY(rot_angle_);
      new G4PVPlacement(G4Transform3D(*window_rot_, window_pos),
                        window_logic_vol, window_name, lab_logic,
                        false, 0, false);


      sensor_x_pos = sensor_x_pos + window_thickness + 1.5 * mm;

    }

    G4ThreeVector sensor_pos = G4ThreeVector(sensor_x_pos, 0., box_z_);

    G4RotationMatrix* sensor_rot_ = new G4RotationMatrix();
    rot_angle_ = -pi/2.;
    // rot_angle_ = 0.;
    sensor_rot_->rotateY(rot_angle_);
    new G4PVPlacement(G4Transform3D(*sensor_rot_, sensor_pos), photo_sensor_logic,
                      photo_sensor_logic->GetName(),lab_logic,true,0,true);


    // Al DISK____________________________________________________

    G4String disk_name = "Al disk";

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


    // BLACK BOX_____________________________________________________

    G4String box_name = "Black box";
    G4String box_inner_name = "inner Black box";

    G4Box* box_outer_solid_vol =
      new G4Box(box_name, box_xy_/2., box_xy_/2., box_z_/2.);

    G4Box* box_inner_solid_vol =
    // new G4Box(box_side_name, box_xy_/2., box_xy_/2., side_thickness_/2.);
    new G4Box(box_inner_name, box_xy_/2. - side_thickness_,
      box_xy_/2. - side_thickness_, box_z_/2. - side_thickness_/2.);

    G4ThreeVector inner_pos = G4ThreeVector(0., 0., side_thickness_/2.);
    G4RotationMatrix* inner_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    rot_angle_ = 0.;
    inner_rot_->rotateY(rot_angle_);
    G4SubtractionSolid* box_solid_vol =
      new G4SubtractionSolid("Box-Side", box_outer_solid_vol, box_inner_solid_vol,
                              inner_rot_, inner_pos);

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

    // TEFLON PANEL_____________________________________________________

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
      panel_opsur->SetMaterialPropertiesTable(opticalprops::PTFE());
    new G4LogicalSkinSurface("Al_OPSURF", panel_logic_vol, panel_opsur);

    G4double panel_z_pos = box_z_ + radius_ + panel_thickness/2.;
    G4ThreeVector panel_pos = G4ThreeVector(0., 0., panel_z_pos);

    G4RotationMatrix* panel_rot_ = new G4RotationMatrix();
    // rot_angle_ = pi;
    rot_angle_ = 0.;
    panel_rot_->rotateY(rot_angle_);
    new G4PVPlacement(G4Transform3D(*panel_rot_, panel_pos),
                      panel_logic_vol, panel_name, lab_logic,
                      false, 0, false);

    // LOOP_____________________________________________________________

    for (int i=0; i < n_fibers; i++){

      // fiber

      x = length_/2 - box_xy_/2 - 1*mm;
      // y = i*radius_ - box_xy_/2;
      y = (box_xy_/n_fibers)*i - box_xy_/2 + radius_/2.;
      z = box_z_ + radius_/2;

      G4RotationMatrix* fib_rot_ = new G4RotationMatrix();
      rot_angle_ = pi/2.;
      // rot_angle_ = 0.;
      fib_rot_->rotateY(rot_angle_);

      new G4PVPlacement(fib_rot_,G4ThreeVector(x, y, z),fiber_logic,
                              fiber_logic->GetName(),lab_logic,true,0,true);


      // Optical gel

      G4ThreeVector opt_gel_pos = G4ThreeVector(x + length_/2 + opt_gel_thickness/2., y, z);
      G4RotationMatrix* opt_gel_rot_ = new G4RotationMatrix();
      rot_angle_ = pi/2.;
      // rot_angle_ = 0.;
      opt_gel_rot_->rotateY(rot_angle_);
      new G4PVPlacement(G4Transform3D(*opt_gel_rot_, opt_gel_pos),
                        opt_gel_logic_vol, opt_gel_name, lab_logic,
                        false, 0, false);


      // Al disk

      G4ThreeVector disk_pos = G4ThreeVector(x - length_/2 - disk_thickness/2., y, z);
      G4RotationMatrix* disk_rot_ = new G4RotationMatrix();
      rot_angle_ = pi/2.;
      // rot_angle_ = 0.;
      disk_rot_->rotateY(rot_angle_);
      new G4PVPlacement(G4Transform3D(*disk_rot_, disk_pos),
                        disk_logic_vol, disk_name, lab_logic,
                        false, 0, false);

    }



}
G4ThreeVector Fib_box_struct::GenerateVertex(const G4String& region) const {
    // return cyl_vertex_gen_->GenerateVertex(region);

    // // G4ThreeVector vertex(1.,1.,1.);
    // G4ThreeVector vertex(box_xy_/2, box_xy_/2, 0.);
    G4ThreeVector vertex(0., 0., side_thickness_);

    // WORLD
    if (region == "CENTER") {
      return vertex;
    }

    return vertex;

}
