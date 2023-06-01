#include "Fiber_barrel_meth.h"
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"
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

REGISTER_CLASS(Fiber_barrel_meth,GeometryBase)

Fiber_barrel_meth::Fiber_barrel_meth():
GeometryBase(),
radius_(1.*mm),
length_(1.*cm),
radius_cyl_(1. *cm),
sensor_type_ ("PERFECT"),
sensor_visibility_ (true),
cyl_vertex_gen_(0)
{
    msg_=new G4GenericMessenger(this,"/Geometry/Fiber_barrel_meth/","Control commands of geometry OpticalFibre");

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

    msg_->DeclareProperty("sensor_type", sensor_type_,
        "Sensors type");

    msg_->DeclareProperty("sensor_visibility", sensor_visibility_,
                          "Sensors visibility");

    cyl_vertex_gen_ = new CylinderPointSampler(0.5*radius_cyl_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);
}
Fiber_barrel_meth::~Fiber_barrel_meth() {
    delete cyl_vertex_gen_;
    delete msg_;
}
void Fiber_barrel_meth::Construct(){
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


    G4double x;
    G4double y;
    G4double z = 0.;
    G4double rot_angle;

    G4Material* ps = materials::PS();
    G4Material* tpb = materials::TPB();


    // Outer teflon cilynder

    G4String teflon_name = "Teflon_panel";

    G4double teflon_thickness = .1 * mm;

    G4Material* teflon_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    teflon_mat->SetMaterialPropertiesTable(opticalprops::PTFE());

    G4Tubs* teflon_solid_vol =
      new G4Tubs(teflon_name, radius_cyl_ + radius_/2., radius_cyl_ + radius_/2. + teflon_thickness, length_/2., 0., 360.*deg);

    G4LogicalVolume* teflon_logic_vol =
      new G4LogicalVolume(teflon_solid_vol, teflon_mat, teflon_name);

    G4OpticalSurface* teflon_opsur =
    new G4OpticalSurface("teflon_OPSURF", unified, polished, dielectric_metal);
    teflon_opsur->SetMaterialPropertiesTable(opticalprops::PTFE());
    new G4LogicalSkinSurface("teflon_OPSURF", teflon_logic_vol, teflon_opsur);

    G4VisAttributes teflon_col = nexus::Lilla();
    teflon_logic_vol->SetVisAttributes(teflon_col);
    // teflon_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4ThreeVector teflon_pos = G4ThreeVector(0., 0., 0.);
    G4RotationMatrix* teflon_rot = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    teflon_rot->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*teflon_rot, teflon_pos),
                      teflon_logic_vol, teflon_name, lab_logic,
                      false, 0, false);



    // Inner methacrylate cilynder

    G4String window_name = "Methacrylate window";

    G4double window_thickness = 1.5 * mm;

    G4Material* window_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
    window_mat->SetMaterialPropertiesTable(opticalprops::PMMA());

    G4Tubs* window_solid_vol =
    new G4Tubs(window_name, radius_cyl_ - radius_/2. - window_thickness, radius_cyl_ - radius_/2., length_/2., 0., 360.*deg);

    G4LogicalVolume* window_logic_vol =
      new G4LogicalVolume(window_solid_vol, window_mat, window_name);
    G4VisAttributes window_col = nexus::LightBlue();
    window_logic_vol->SetVisAttributes(window_col);

    G4OpticalSurface* window_opsur =
      new G4OpticalSurface("window_OPSURF", unified, polished, dielectric_dielectric);
      // window_opsur->SetMaterialPropertiesTable(opticalprops::PMMA());
      window_opsur->SetMaterialPropertiesTable(opticalprops::TPB());
    new G4LogicalSkinSurface("window_OPSURF", window_logic_vol, window_opsur);

    G4ThreeVector window_pos = G4ThreeVector(0., 0., 0.);

    G4RotationMatrix* window_rot_ = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    window_rot_->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*window_rot_, window_pos),
                      window_logic_vol, window_name, lab_logic,
                      false, 0, false);

    // Al disk to aluminize fibers _____________________________________________________________________

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


   // Sensors at the not aluminized end of the fibers _____________________________________________________________________

   G4double sensor_thickness = 1. * mm;

    /// Build the sensor
    photo_sensor_  = new GenericPhotosensor("F_SENSOR", radius_,
                                            radius_, sensor_thickness);
    /// Constructing the sensors
    // Optical Properties of the sensor
    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();

    if (sensor_type_ == "PERFECT") {
      // perfect detector
      G4int entries = 4;
      G4double energy[entries]       = {0.2 * eV, 3.5 * eV, 3.6 * eV, 11.5 * eV};
      G4double efficiency[entries]   = {1.      , 1.      , 1.      , 1.       };

      photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   entries);
    }

    else if (sensor_type_ == "PMT") {
      // PMT
      G4int entries = 9;
      G4double energy[entries]       = {
        h_Planck * c_light / (903.715 * nm), h_Planck * c_light / (895.975 * nm),
        h_Planck * c_light / (866.563 * nm), h_Planck * c_light / (826.316 * nm),
        h_Planck * c_light / (628.173 * nm), h_Planck * c_light / (490.402 * nm),
        h_Planck * c_light / (389.783 * nm), h_Planck * c_light / (330.96 * nm),
        h_Planck * c_light / (296.904 * nm)
      };
      G4double efficiency[entries]   = {0.00041, 0.00107, 0.01248, 0.06181,
                                        0.12887, 0.19246, 0.09477, 0.06040,
                                        0.00826};

      photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   entries);
    }

    else if (sensor_type_ == "SiPM_FBK") {
      // SiPM_FBK
      G4int entries = 13;
      G4double energy[entries]       = {
        h_Planck * c_light / (699.57 * nm), h_Planck * c_light / (630.00 * nm),
        h_Planck * c_light / (590.43 * nm), h_Planck * c_light / (544.78 * nm),
        h_Planck * c_light / (524.78 * nm), h_Planck * c_light / (499.57 * nm),
        h_Planck * c_light / (449.57 * nm), h_Planck * c_light / (435.22 * nm),
        h_Planck * c_light / (420.00 * nm), h_Planck * c_light / (409.57 * nm),
        h_Planck * c_light / (399.57 * nm), h_Planck * c_light / (389.57 * nm),
        h_Planck * c_light / (364.78 * nm)
      };
      G4double efficiency[entries]   = {.2137, .2743,
                                        .3189, .3634,
                                        .3829, .4434,
                                        .4971, .5440,
                                        .5657, .5829,
                                        .5886, .5657,
                                        .4743
                                        };

      photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   entries);
    }

    else if (sensor_type_ == "SiPM_Hamamatsu") {
      // SiPM_Hamamatsu
      G4int entries = 13;
      G4double energy[entries]       = {
        h_Planck * c_light / (831.818 * nm), h_Planck * c_light / (761.932 * nm),
        h_Planck * c_light / (681.818 * nm), h_Planck * c_light / (620.455 * nm),
        h_Planck * c_light / (572.727 * nm), h_Planck * c_light / (516.477 * nm),
        h_Planck * c_light / (460.227 * nm), h_Planck * c_light / (400.568 * nm),
        h_Planck * c_light / (357.955 * nm), h_Planck * c_light / (344.318 * nm),
        h_Planck * c_light / (311.932 * nm), h_Planck * c_light / (289.773 * nm),
        h_Planck * c_light / (282.955 * nm)
      };
      G4double efficiency[entries]   = {.07329, .12673,
                                        .20254, .29851,
                                        .36889, .45739,
                                        .49695, .44929,
                                        .35476, .35374,
                                        .29960, .19862,
                                        .12204
                                        };

      photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   entries);
    }


    G4double MinE_MaxE[] = {0.2 * eV, 11.5 * eV};
    G4double reflectivity[] = {0., 0.};
    photosensor_mpt->AddProperty("REFLECTIVITY", MinE_MaxE, reflectivity, 2);

    photo_sensor_ ->SetOpticalProperties(photosensor_mpt);

    // Adding to sensors encasing, the Refractive Index of fibers to avoid reflections
    std::cout<<"HERE!"<<std::endl;

    G4Material* fiber_mat = materials::Y11();
    fiber_mat->SetMaterialPropertiesTable(opticalprops::Y11());

    G4MaterialPropertyVector* fibers_rindex =
    fiber_mat->GetMaterialPropertiesTable()->GetProperty("RINDEX");

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


    // Fibers _____________________________________________________________________

    G4double n_fibers = radius_cyl_*pi/radius_;
    G4double dif_theta = .5*2*pi/n_fibers; // angular separation between fibers
    G4double theta;


    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;

    // fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true); // coated
    fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, false, tpb, ps, true); // uncoated
    fiber_->SetCoreOpticalProperties(opticalprops::Y11());
    // fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
    fiber_->Construct();
    fiber_logic = fiber_->GetLogicalVolume();


    // Loop to place elements

    for (int i=0; i<=2*n_fibers; i++){

      theta = dif_theta*i;
      x = radius_cyl_ * cos(theta);
      y = radius_cyl_ * sin(theta);

      // fibers
      new G4PVPlacement(0,G4ThreeVector(x, y, z),fiber_logic,
                        fiber_logic->GetName(),lab_logic,true,0,true);

      // aluminization
      new G4PVPlacement(0,G4ThreeVector(x, y, z - length_/2. - disk_thickness/2.),
                        disk_logic_vol, disk_name,lab_logic,true,0,true);

      // sensor
      G4ThreeVector sensor_pos = G4ThreeVector(x, y, z + length_/2. + sensor_thickness/2.);

      G4RotationMatrix* sensor_rot = new G4RotationMatrix();
      // rot_angle = -pi/2.;
      // rot_angle = 0.;
      rot_angle = pi;
      sensor_rot->rotateY(rot_angle);
      new G4PVPlacement(G4Transform3D(*sensor_rot, sensor_pos), photo_sensor_logic,
                        photo_sensor_logic->GetName(),lab_logic,true,0,true);



    }


}


G4ThreeVector Fiber_barrel_meth::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);

    // // G4ThreeVector vertex(1.,1.,1.);
    // G4ThreeVector vertex(0., 0., 0.);
    // return vertex;

}
