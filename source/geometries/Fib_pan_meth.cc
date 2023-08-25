// ----------------------------------------------------------------------------
// nexus | Fib_pan_meth.cc
//
// Box containing optical fibers
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Fib_pan_meth.h"

#include "FactoryBase.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
#include "IonizationSD.h"

#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4UserLimits.hh>

using namespace nexus;

REGISTER_CLASS(Fib_pan_meth, GeometryBase)

namespace nexus {

  Fib_pan_meth::Fib_pan_meth():
    GeometryBase(),
    liquid_(false),
    pressure_(STP_Pressure),
    radius_ (1 * cm),
    methacrylate_ (false),
    window_thickness_ (5. * mm),
    fiber_diameter_(1 * mm),
    length_ (2 * cm),
    fiber_type_ ("Y11"),
    sensor_type_ ("PERFECT"),
    sensor_visibility_ (true),
    teflon_thickness_ (1. *mm),
    panel_width_ (170. * mm),
    caps_visibility_ (false),
    teflon_visibility_ (false),
    coated_(true)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/Fib_pan_meth/",
      "Control commands of geometry Fib_pan_meth.");

    msg_->DeclareProperty("LXe", liquid_,
      "Fill the inside of the detector with liquid xenon.");

    G4GenericMessenger::Command& pressure_cmd =
      msg_->DeclareProperty("pressure", pressure_,
      "Set pressure for gaseous xenon (if selected).");
    pressure_cmd.SetUnitCategory("Pressure");
    pressure_cmd.SetParameterName("pressure", false);
    pressure_cmd.SetRange("pressure>0.");

    G4GenericMessenger::Command&  radius_cmd =
      msg_->DeclareProperty("radius", radius_,
                            "Barrel fiber radius");
    radius_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  window_thickness_cmd =
      msg_->DeclareProperty("window_thickness", window_thickness_,
                            "Inner methacrylate thickness");
    window_thickness_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  length_cmd =
      msg_->DeclareProperty("length", length_,
                            "Barrel fiber length");
    length_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  fiber_diameter_cmd =
      msg_->DeclareProperty("fiber_diameter", fiber_diameter_,
                            "Fiber diameter");
    fiber_diameter_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  teflon_thickness_cmd =
      msg_->DeclareProperty("teflon_thickness", teflon_thickness_,
                            "Teflon coat thickness");
    teflon_thickness_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  panel_width_cmd =
      msg_->DeclareProperty("panel_width", panel_width_,
                            "Teflon panel width");
    panel_width_cmd.SetUnitCategory("Length");

    msg_->DeclareProperty("methacrylate", methacrylate_, "Inner methacrylate window (true or false)");
    msg_->DeclareProperty("fiber_type", fiber_type_, "Fiber type (Y11 or B2)");
    msg_->DeclareProperty("sensor_type", sensor_type_, "Sensors type");
    msg_->DeclareProperty("sensor_visibility", sensor_visibility_, "Sensors visibility");
    msg_->DeclareProperty("caps_visibility", caps_visibility_, "Make teflon endcaps visible (true or false)");
    msg_->DeclareProperty("teflon_visibility", teflon_visibility_, "Make teflon coat visible (true or false)");
    msg_->DeclareProperty("coated", coated_, "Coat fibers with WLS coating");

  }



  Fib_pan_meth::~Fib_pan_meth()
  {
    delete msg_;
  }



  void Fib_pan_meth::Construct()
  {

    G4cout << "[Fib_pan_meth] *** Barrel Fiber prototype ***" << G4endl;
    G4cout << "[Fib_pan_meth] Using " << fiber_type_ << " fibers";


    world_z_ = (length_ + teflon_thickness_) * 3;
    world_xy_ = (radius_ + teflon_thickness_) * 3;

    G4double rot_angle;


    G4Material *this_fiber = nullptr;
    G4MaterialPropertiesTable *this_fiber_optical = nullptr;
    G4Material *this_coating = nullptr;
    G4MaterialPropertiesTable *this_coating_optical = nullptr;

    if (fiber_type_ == "Y11") {
      this_fiber = materials::Y11();
      this_fiber_optical = opticalprops::Y11();

      if (coated_) {
        this_coating = materials::TPB();
        this_coating_optical = opticalprops::TPB();
      }

    } else if (fiber_type_ == "B2") {

      this_fiber = materials::B2();
      this_fiber_optical = opticalprops::B2();

      if (coated_) {
        this_coating = materials::TPH();
        this_coating_optical = opticalprops::TPH();
      }

    } else {
      G4Exception("[FiberBarrel]", "Construct()",
                  FatalException, "Invalid fiber type, must be Y11 or B2");
    }


    fiber_ = new GenericWLSFiber(fiber_type_, true, fiber_diameter_, length_, true, coated_, this_coating, this_fiber, true);

    // WORLD /////////////////////////////////////////////////

    // Define the material (LXe or GXe) for the tank.
    // We use for this the NIST manager or the nexus materials list.
    G4Material* xenon = 0;
    if (liquid_)
      xenon = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");
    else
      xenon = materials::GXe(pressure_);

    G4String world_name = "WORLD";

    // G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* world_mat = xenon;

    // world_mat->SetMaterialPropertiesTable(opticalprops::Vacuum());
    world_mat->SetMaterialPropertiesTable(opticalprops::GXe());

    G4Box* world_solid_vol =
     new G4Box(world_name, world_xy_/2., world_xy_/2., world_z_/2.);

    G4LogicalVolume* world_logic_vol =
      new G4LogicalVolume(world_solid_vol, world_mat, world_name);
    world_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());
    GeometryBase::SetLogicalVolume(world_logic_vol);

    // // Set the logical volume of the lab as an ionization
    // // sensitive detector, i.e. position, time and energy deposition
    // // will be stored for each step of any charged particle crossing
    // // the volume.
    // IonizationSD* ionizsd = new IonizationSD("/Fib_pan_meth");
    // G4SDManager::GetSDMpointer()->AddNewDetector(ionizsd);
    // world_logic_vol->SetSensitiveDetector(ionizsd);
    // // Fix the length of the maximum step an electron can make
    // // when depositing energy.
    // // The smaller the limit the narrower the trace
    // world_logic_vol->SetUserLimits(new G4UserLimits(1.*mm));



    // GEOMETRY PARAMETERS /////////////////////////////////////////////

    // Teflon panels distance to the center
    // G4double h = radius_ + fiber_diameter_/2. + teflon_thickness_/2.;
    G4double h = radius_ - teflon_thickness_/2.;

    // Teflon panels angular separation
    G4double dif_theta = 2*std::atan(panel_width_/(2.*h));

    G4int n_panels = floor(( 2 * M_PI) / dif_theta); // optimize the number of panels

    G4int n_fibers = floor(panel_width_ / fiber_diameter_); // number of fibers per panel
    G4double dl_fib = panel_width_/n_fibers; // distance between fibers

    G4int n_sensors = 5; // number of sensors per panel
    G4double dl_sens = panel_width_/n_sensors; // distance between sensors

    // Re-calculation of parameters for optimization
    dif_theta = ( 2 * M_PI) / n_panels; // re-calculate angular difference
    h = (panel_width_/2.)/(std::tan(dif_theta/2.)) + teflon_thickness_/2. + fiber_diameter_; // re-calculate distance to the center
    G4cout << "[Fib_pan_meth] Using " << n_panels << " panels" << G4endl;

    // Fibers/sensors/aluminium distance to the center
    G4double hh = h - (fiber_diameter_/2. + teflon_thickness_/2.);


    // volume where the vertex will be created
    if (methacrylate_)
    {
      cyl_vertex_gen_ = new CylinderPointSampler2020(0, hh - fiber_diameter_/2. - window_thickness_, length_/2, 0, 2 * M_PI);
    } else if (methacrylate_ == false){
      cyl_vertex_gen_ = new CylinderPointSampler2020(0, hh - fiber_diameter_/2., length_/2, 0, 2 * M_PI,
                                                     nullptr, G4ThreeVector(0,0,0));
      // cyl_vertex_gen_ = new CylinderPointSampler2020(0, hh/2., length_/2, 0, 2 * M_PI);
    }


    // TEFLON END-CAPS ////////////////////////////////////////////
    G4Tubs* teflon_cap =
    // new G4Tubs("TEFLON_CAP", 0, radius_ - (fiber_diameter_ + teflon_thickness_) , teflon_thickness_/2., 0, twopi);
    new G4Tubs("TEFLON_CAP", 0, hh - fiber_diameter_/2., teflon_thickness_/2., 0, twopi);
    G4Material* teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    teflon->SetMaterialPropertiesTable(opticalprops::PTFE());
    G4LogicalVolume* teflon_logic =
      new G4LogicalVolume(teflon_cap, teflon, "TEFLON");

    G4OpticalSurface* opsur_teflon =
      new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
    opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());

    new G4LogicalSkinSurface("TEFLON_OPSURF", teflon_logic, opsur_teflon);

    if (caps_visibility_ == false)
    {
      teflon_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, -(length_ + teflon_thickness_)/2.),
                      teflon_logic, "TEFLON1", world_logic_vol,
                      true, 0, false);

    new G4PVPlacement(0, G4ThreeVector(0, 0, (length_ + teflon_thickness_)/2.),
                      teflon_logic, "TEFLON2", world_logic_vol,
                      true, 1, false);

    // TEFLON PANEL /////////////////////////////////////////////////

    G4Box* teflon_panel =
      new G4Box("TEFLON_PANEL", panel_width_/2., length_/2., teflon_thickness_/2.);
    G4LogicalVolume* teflon_panel_logic =
      new G4LogicalVolume(teflon_panel, teflon, "TEFLON_PANEL");

    new G4LogicalSkinSurface("TEFLON_PANEL_OPSURF", teflon_panel_logic, opsur_teflon);

    teflon_panel_logic->SetVisAttributes(nexus::White());

    if (teflon_visibility_ == false)
    {
      teflon_panel_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }

    // FIBER ////////////////////////////////////////////////////

    fiber_->SetCoreOpticalProperties(this_fiber_optical);
    fiber_->SetCoatingOpticalProperties(this_coating_optical);

    fiber_->Construct();
    G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
    if (fiber_type_ == "Y11")
      fiber_logic->SetVisAttributes(nexus::LightGreenAlpha());
    else if (fiber_type_ == "B2")
      fiber_logic->SetVisAttributes(nexus::LightBlueAlpha());

    // G4int n_fibers = floor((radius_ * 2 * M_PI) / fiber_diameter_);
    // G4cout << "[Fib_pan_meth] Barrel with " << n_fibers << " fibers" << G4endl;

    // DETECTOR /////////////////////////////////////////////////
    G4double sensor_width = panel_width_/n_sensors;
    G4double sensor_thickness = 1. * mm;
    G4String sensor_name = "F_SENSOR";

    /// Build the sensor
    photo_sensor_  = new GenericPhotosensor(sensor_name, sensor_width, fiber_diameter_, sensor_thickness);


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


      G4Material* window_mat = this_fiber;
      window_mat->SetMaterialPropertiesTable(this_fiber_optical);

      G4MaterialPropertyVector* window_rindex =
      window_mat->GetMaterialPropertiesTable()->GetProperty("RINDEX");


      photo_sensor_ ->SetOpticalProperties(photosensor_mpt);

      // Adding to sensors encasing, the Refractive Index of fibers to avoid reflections
      photo_sensor_ ->SetWindowRefractiveIndex(window_rindex);

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



    // ALUMINIZED ENDCAP//////////////////////////////////////////////////

    G4Material* fiber_end_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    G4double fiber_end_z = 0.1 * mm;

    G4Tubs* fiber_end_solid_vol =
      new G4Tubs("fiber_end", 0, fiber_diameter_ / 2, fiber_end_z/2., 0, 2 * M_PI);

    G4LogicalVolume* fiber_end_logic_vol =
      new G4LogicalVolume(fiber_end_solid_vol, fiber_end_mat, "FIBER_END");
    G4OpticalSurface* opsur_al =
      new G4OpticalSurface("POLISHED_AL_OPSURF", unified, polished, dielectric_metal);
    opsur_al->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("POLISHED_AL_OPSURF", fiber_end_logic_vol, opsur_al);

    fiber_end_logic_vol  ->SetVisAttributes(nexus::Blue());


    // PLACEMENT /////////////////////////////////////////////

    for (G4int itheta=0; itheta < n_panels; itheta++) {
    // for (G4int itheta=0; itheta < 3; itheta++) {

      // panels
      G4double theta = dif_theta * itheta;
      G4double x = h * std::cos(theta) * mm;
      G4double y = h * std::sin(theta) * mm;
      G4double phi = pi/2. + std::atan2(y, x);
      std::string label = std::to_string(itheta);

      G4RotationMatrix* panel_rot = new G4RotationMatrix();
      rot_angle = pi/2.;
      panel_rot->rotateX(rot_angle);
      panel_rot->rotateY(phi);
      new G4PVPlacement(panel_rot, G4ThreeVector(x,y, 0.),
                        teflon_panel_logic, "PANEL-"+label, world_logic_vol,
                        false, itheta, false);

      // INNER METHACRYLATE PANELS/////////////////////////////////////////////
      if (methacrylate_) {

        G4String window_name = "METHACRYLATE_WINDOW";

        G4Material* window_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
        window_mat->SetMaterialPropertiesTable(opticalprops::PMMA());

        G4Box* window_solid_vol =
          new G4Box(window_name, panel_width_/2., length_/2., window_thickness_/2.);

        G4LogicalVolume* window_logic_vol =
          new G4LogicalVolume(window_solid_vol, window_mat, window_name);
        G4VisAttributes window_col = nexus::Lilla();
        window_logic_vol->SetVisAttributes(window_col);

        G4OpticalSurface* window_opsur =
          new G4OpticalSurface("window_OPSURF", unified, polished, dielectric_dielectric);
          // window_opsur->SetMaterialPropertiesTable(opticalprops::PMMA());
          if (fiber_type_ == "Y11") window_opsur->SetMaterialPropertiesTable(opticalprops::TPB());
          if (fiber_type_ == "B2") window_opsur->SetMaterialPropertiesTable(opticalprops::TPH());
        new G4LogicalSkinSurface("window_OPSURF", window_logic_vol, window_opsur);

        G4double hw = h - (window_thickness_/2. + fiber_diameter_ + teflon_thickness_/2.);

        new G4PVPlacement(panel_rot, G4ThreeVector(x*hw/h,y*hw/h, 0.),
                          window_logic_vol, window_name + label, world_logic_vol,
                          false, 0, false);
      }

      // Relative positions of the fibers wrt the panel
      G4double x0_f = x*hh/h + (panel_width_/2. - fiber_diameter_/2.)*std::cos(phi);
      G4double y0_f = y*hh/h + (panel_width_/2. - fiber_diameter_/2.)*std::sin(phi);

      // Relative positions of the sensors wrt the panel
      G4double x0_s = x*hh/h + (panel_width_/2. - sensor_width/2.)*std::cos(phi);
      G4double y0_s = y*hh/h + (panel_width_/2. - sensor_width/2.)*std::sin(phi);


      for (G4int ii=0; ii < n_fibers; ii++) {
      // for (G4int ii=0; ii < 1; ii++) {

          G4double xx_f = x0_f - dl_fib*ii*std::cos(phi);
          G4double yy_f = y0_f - dl_fib*ii*std::sin(phi);

          std::string label2 = std::to_string(ii);

          new G4PVPlacement(0, G4ThreeVector(xx_f, yy_f),
                            fiber_logic, "FIBER-"+label+label2, world_logic_vol,
                            false, n_panels + ii, false);
          new G4PVPlacement(0, G4ThreeVector(xx_f, yy_f, -(length_ + fiber_end_z)/2.),
                            fiber_end_logic_vol, "ALUMINIUMR-" + label + label2, world_logic_vol,
                            false, 2*n_panels + ii, false);
        }
      for (G4int jj=0; jj < n_sensors; jj++) {
      // for (G4int jj=0; jj < 3; jj++) {

            G4double xx_s = x0_s - dl_sens*jj*std::cos(phi);
            G4double yy_s = y0_s - dl_sens*jj*std::sin(phi);

            std::string label2 = std::to_string(jj);

            G4RotationMatrix* sensor_rot = new G4RotationMatrix();
            // rot_angle = 0.;
            rot_angle = M_PI;
            sensor_rot->rotateY(rot_angle);
            sensor_rot->rotateZ(phi);
            new G4PVPlacement(sensor_rot, G4ThreeVector(xx_s, yy_s, (length_ + sensor_thickness)/2.),
                              photo_sensor_logic, "SENS-" + label+label2, world_logic_vol,
                              true, 3*n_panels + jj, false);

      }

    }

  }



  G4ThreeVector Fib_pan_meth::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex;

    if (region == "VOLUME"){

      vertex = cyl_vertex_gen_->GenerateVertex("VOLUME");

    }


    return vertex;
  }


} // end namespace nexus
