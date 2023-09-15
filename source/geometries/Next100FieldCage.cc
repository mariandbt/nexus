// ----------------------------------------------------------------------------
// nexus | Next100FieldCage.cc
//
// Geometry of the NEXT-100 field cage.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100FieldCage.h"
#include "MaterialsList.h"
#include "Visibilities.h"
#include "IonizationSD.h"
#include "OpticalMaterialProperties.h"
#include "UniformElectricDriftField.h"
#include "XenonProperties.h"
#include "CylinderPointSampler2020.h"

#include <G4Navigator.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>
#include <G4GenericMessenger.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Tubs.hh>
#include <G4Polyhedra.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4Box.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4NistManager.hh>
#include <G4UserLimits.hh>
#include <G4SDManager.hh>
#include <G4UnitsTable.hh>
#include <G4TransportationManager.hh>

// Marian's adenda
//

using namespace nexus;


Next100FieldCage::Next100FieldCage():
  GeometryBase(),
  // Dimensions
  active_diam_         (984. * mm), // distance between the centers of two opposite panels

  cathode_int_diam_    (960. * mm),
  cathode_ext_diam_    (1020.* mm),
  cathode_thickn_      (10.  * mm),
  // Caution: updating grid-thickn_ will require updating gate-tp and gate-sapphire-window distances
  grid_thickn_         (0.1  * mm),

  teflon_drift_length_ (1178.*mm), //distance from the gate to the beginning of the cathode volume.
  teflon_total_length_ (1431. * mm),
  teflon_thickn_       (5. * mm),
  // n_panels_            (18),

  el_gap_length_ (10. * mm),

  gate_teflon_dist_ (10.2 * mm - grid_thickn_), //distance from gate-grid to teflon
  gate_ext_diam_    (1042. * mm), //preliminary
  gate_int_diam_    (1009. * mm), //preliminary
  gate_ring_thickn_ (9.9   * mm), // maximum possible value to avoid overlap with sipm board masks

  // external to teflon (hdpe + rings + holders)
  hdpe_tube_int_diam_ (1080. * mm),
  hdpe_tube_ext_diam_ (1105.4 * mm),
  hdpe_length_        (1192. * mm),

  ring_ext_diam_ (1038. * mm),
  ring_int_diam_ (1014. * mm),
  ring_thickn_   (10. * mm),
  drift_ring_dist_  (24. * mm),
  buffer_ring_dist_ (48. * mm),
  holder_x_         (60. * mm),  //x dimension of the holders
  holder_long_y_    (9.  * mm),  // y dim of the base of the ring holders
  holder_short_y_   (33.15 * mm),// y dim of the pieces added over the base of the ring holders

  tpb_thickn_ (1 * micrometer),
  overlap_    (0.001*mm), //defined for G4UnionSolids to ensure a common volume within the two joined solids
  // Diffusion constants
  drift_transv_diff_ (1. * mm/sqrt(cm)),
  drift_long_diff_ (.3 * mm/sqrt(cm)),
  ELtransv_diff_ (0. * mm/sqrt(cm)),
  ELlong_diff_ (0. * mm/sqrt(cm)),
  // EL electric field
  elfield_ (0),
  ELelectric_field_ (34.5*kilovolt/cm),
  cath_grid_transparency_(.95),
  el_grid_transparency_  (.90),
  max_step_size_ (1. * mm),
  visibility_ (0),
  verbosity_(0),
  // EL gap generation disk parameters
  el_gap_gen_disk_diam_(0.),
  el_gap_gen_disk_x_(0.), el_gap_gen_disk_y_(0.),
  el_gap_gen_disk_zmin_(0.), el_gap_gen_disk_zmax_(1.),
  // Fiber Barrel
  fiber_type_ ("Y11"), // type of fibers attached to the teflon panels (Y11 or B2)
  sensor_type_ ("PERFECT"),
  fiber_diameter_(1 * mm),
  panel_width_ (170. * mm),
  sensor_visibility_ (true),
  cap_visibility_ (false),
  panels_visibility_ (false),
  coated_(true)
{
  /// Define new categories
  new G4UnitDefinition("kilovolt/cm","kV/cm","Electric field", kilovolt/cm);
  new G4UnitDefinition("mm/sqrt(cm)","mm/sqrt(cm)","Diffusion", mm/sqrt(cm));

  /// Initializing the geometry navigator (used in vertex generation)
  geom_navigator_ =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  /// Messenger
  msg_ = new G4GenericMessenger(this, "/Geometry/Next100/",
                                "Control commands of geometry Next100.");
  msg_->DeclareProperty("field_cage_vis", visibility_, "Field Cage Visibility");
  msg_->DeclareProperty("field_cage_verbosity", verbosity_, "Field Cage Verbosity");

  G4GenericMessenger::Command& drift_transv_diff_cmd =
    msg_->DeclareProperty("drift_transv_diff", drift_transv_diff_,
                          "Tranvsersal diffusion in the drift region");
  drift_transv_diff_cmd.SetParameterName("drift_transv_diff", true);
  drift_transv_diff_cmd.SetUnitCategory("Diffusion");

  G4GenericMessenger::Command& drift_long_diff_cmd =
  msg_->DeclareProperty("drift_long_diff", drift_long_diff_,
                        "Longitudinal diffusion in the drift region");
  drift_long_diff_cmd.SetParameterName("drift_long_diff", true);
  drift_long_diff_cmd.SetUnitCategory("Diffusion");

  G4GenericMessenger::Command&  ELtransv_diff_cmd =
  msg_->DeclareProperty("ELtransv_diff", ELtransv_diff_,
                        "Tranvsersal diffusion in the EL region");
  ELtransv_diff_cmd.SetParameterName("ELtransv_diff", true);
  ELtransv_diff_cmd.SetUnitCategory("Diffusion");

  G4GenericMessenger::Command&  ELlong_diff_cmd =
  msg_->DeclareProperty("ELlong_diff", ELlong_diff_,
                        "Longitudinal diffusion in the EL region");
  ELlong_diff_cmd.SetParameterName("ELlong_diff", true);
  ELlong_diff_cmd.SetUnitCategory("Diffusion");

  msg_->DeclareProperty("elfield", elfield_,
                        "True if the EL field is on (full simulation), false if it's not (parametrized simulation.");

  G4GenericMessenger::Command& El_field_cmd =
  msg_->DeclareProperty("EL_field", ELelectric_field_,
                        "Electric field in the EL region");
  El_field_cmd.SetParameterName("EL_field", true);
  El_field_cmd.SetUnitCategory("Electric field");

  G4GenericMessenger::Command& step_cmd =
    msg_->DeclareProperty("max_step_size", max_step_size_, "Maximum Step Size");
  step_cmd.SetUnitCategory("Length");
  step_cmd.SetParameterName("max_step_size", true);
  step_cmd.SetRange("max_step_size>0.");

  G4GenericMessenger::Command& el_gap_gen_disk_diam_cmd =
    msg_->DeclareProperty("el_gap_gen_disk_diam", el_gap_gen_disk_diam_,
                          "Diameter of the EL gap vertex generation disk.");
  el_gap_gen_disk_diam_cmd.SetUnitCategory("Length");
  el_gap_gen_disk_diam_cmd.SetParameterName("el_gap_gen_disk_diam", false);
  el_gap_gen_disk_diam_cmd.SetRange("el_gap_gen_disk_diam>=0.");

  G4GenericMessenger::Command& el_gap_gen_disk_x_cmd =
    msg_->DeclareProperty("el_gap_gen_disk_x", el_gap_gen_disk_x_,
                          "X position of the center of the EL gap vertex generation disk.");
  el_gap_gen_disk_x_cmd.SetUnitCategory("Length");
  el_gap_gen_disk_x_cmd.SetParameterName("el_gap_gen_disk_x", false);

  G4GenericMessenger::Command& el_gap_gen_disk_y_cmd =
    msg_->DeclareProperty("el_gap_gen_disk_y", el_gap_gen_disk_y_,
                          "Y position of the center of the EL gap vertex generation disk.");
  el_gap_gen_disk_y_cmd.SetUnitCategory("Length");
  el_gap_gen_disk_y_cmd.SetParameterName("el_gap_gen_disk_y", false);

  G4GenericMessenger::Command& el_gap_gen_disk_zmin_cmd =
    msg_->DeclareProperty("el_gap_gen_disk_zmin", el_gap_gen_disk_zmin_,
                          "Minimum Z range of the EL gap vertex generation disk.");
  el_gap_gen_disk_zmin_cmd.SetParameterName("el_gap_gen_disk_zmin", false);
  el_gap_gen_disk_zmin_cmd.SetRange("el_gap_gen_disk_zmin>=0.0 && el_gap_gen_disk_zmin<=1.0");

  G4GenericMessenger::Command& el_gap_gen_disk_zmax_cmd =
    msg_->DeclareProperty("el_gap_gen_disk_zmax", el_gap_gen_disk_zmax_,
                          "Maximum Z range of the EL gap vertex generation disk.");
  el_gap_gen_disk_zmax_cmd.SetParameterName("el_gap_gen_disk_zmax", false);
  el_gap_gen_disk_zmax_cmd.SetRange("el_gap_gen_disk_zmax>=0.0 && el_gap_gen_disk_zmax<=1.0");

  // Fiber Barrel
  G4GenericMessenger::Command&  fiber_diameter_cmd =
      msg_->DeclareProperty("fiber_diameter", fiber_diameter_,
                            "Fiber diameter");
  fiber_diameter_cmd.SetUnitCategory("Length");

  G4GenericMessenger::Command&  panel_width_cmd =
    msg_->DeclareProperty("panel_width", panel_width_,
                          "Teflon panel width");
  panel_width_cmd.SetUnitCategory("Length");

  msg_->DeclareProperty("fiber_type", fiber_type_, "Fiber type (Y11 or B2)");
  msg_->DeclareProperty("sensor_type", sensor_type_, "Sensors type");
  msg_->DeclareProperty("sensor_visibility", sensor_visibility_, "Sensors visibility");
  msg_->DeclareProperty("cap_visibility", cap_visibility_, "Make teflon endcap visible (true or false)");
  msg_->DeclareProperty("panels_visibility", panels_visibility_, "Make teflon panels visible (true or false)");
  msg_->DeclareProperty("coated", coated_, "Coat fibers with WLS coating");


}


void Next100FieldCage::SetMotherLogicalVolume(G4LogicalVolume* mother_logic)
{
  mother_logic_ = mother_logic;
}


void Next100FieldCage::SetMotherPhysicalVolume(G4VPhysicalVolume* mother_phys)
{
  mother_phys_ = mother_phys;
}


void Next100FieldCage::Construct()
{
  // FIBER BARREL ******************************************************
  /// DIMENSIONS PARAMETERS /////////////////////////////////////////////
  panel_thickness_ = teflon_thickn_;
  panel_length_ = teflon_drift_length_;

  sens_z = 1. * mm;
  fiber_end_z = 0.1 * mm;
  fiber_length = panel_length_ - (sens_z + fiber_end_z);
  // ****************************************************************************

  /// Calculate lengths of active and buffer regions
  active_length_ = (cathode_thickn_ - grid_thickn_)/2. + teflon_drift_length_ + gate_teflon_dist_;
  buffer_length_ = gate_sapphire_wdw_dist_ - active_length_ - grid_thickn_;

  /// Calculate length of teflon in the buffer region
  teflon_buffer_length_ = teflon_total_length_ - cathode_thickn_ - teflon_drift_length_;

  /// Calculate radial position of the ring holders.
  holder_r_ = (active_diam_+2.*teflon_thickn_+holder_long_y_)/2.;

  /// Calculate relative positions in mother volume
  gate_grid_zpos_  = GetELzCoord() - grid_thickn_/2.;
  active_zpos_     = gate_grid_zpos_ + grid_thickn_/2. + active_length_/2.;
  cathode_zpos_    = gate_grid_zpos_ + grid_thickn_/2. + active_length_ + grid_thickn_/2.;
  gate_zpos_       = gate_grid_zpos_ + grid_thickn_/2. + gate_ring_thickn_/2. - grid_thickn_;
  el_gap_zpos_     = gate_grid_zpos_ - grid_thickn_/2. - el_gap_length_/2.;
  anode_zpos_      = el_gap_zpos_ - el_gap_length_/2. - gate_ring_thickn_/2.;
  anode_grid_zpos_ = anode_zpos_ + gate_ring_thickn_/2. - grid_thickn_/2.;

  teflon_drift_zpos_  = gate_grid_zpos_ + grid_thickn_/2. + gate_teflon_dist_ + teflon_drift_length_/2.;
  teflon_buffer_zpos_ = cathode_zpos_ + cathode_thickn_/2. + teflon_buffer_length_/2.;

  if (verbosity_) {
    G4cout << "Active length = " << active_length_/mm << " mm" << G4endl;
    G4cout << "Buffer length = " << buffer_length_/mm << " mm" << G4endl;
    G4cout << G4endl;
  }

  // FIBER BARREL ******************************************************
  /// GEOMETRY PARAMETERS /////////////////////////////////////////////

  z = gate_grid_zpos_ + grid_thickn_/2. + active_length_/2.; // z-position of the panels
  z_f = z + fiber_length/2. + sens_z - panel_length_/2.; // z-pos for the fibers
  z_fend = z_f + (fiber_length + fiber_end_z)/2.; // z-pos for the fibers' Al ends
  z_s = z_f - (fiber_length + sens_z)/2.; // z-pos for the sensors

  //// Teflon panels distance to the center
  h = active_diam_/2. - panel_thickness_/2.;

  //// Teflon panels angular separation
  dif_theta = 2*std::atan(panel_width_/(2.*h));

  n_panels = floor(( 2 * M_PI) / dif_theta); // optimize the number of panels

  n_fibers = floor(panel_width_ / fiber_diameter_); // number of fibers per panel
  dl_fib = panel_width_/n_fibers; // distance between fibers

  G4cout << "[FiberBarrel] Using " << n_fibers << " fibers per panel"<< G4endl;


  n_sensors = 5; // number of sensors per panel
  dl_sens = panel_width_/n_sensors; // distance between sensors

  G4cout << "[FiberBarrel] Using " << n_panels*n_sensors << " sensors in total"<< G4endl;


  //// Re-calculation of parameters for optimization
  dif_theta = ( 2 * M_PI) / n_panels; // re-calculate angular difference
  h = (panel_width_/2.)/(std::tan(dif_theta/2.)) + panel_thickness_/2. + fiber_diameter_; // re-calculate distance to the center
  G4cout << "[FiberBarrel] Using " << n_panels << " panels" << G4endl;

  //// Fibers/sensors/aluminium distance to the center
  hh = h - (fiber_diameter_/2. + panel_thickness_/2.);

 // ****************************************************************************


  /// Define materials to be used
  DefineMaterials();
  /// Build the different parts of the field cage
  BuildActive();
  BuildCathode();
  BuildBuffer();
  BuildELRegion();
  BuildFiberBarrel();
  // BuildLightTube();
  BuildFieldCage();
}


void Next100FieldCage::DefineMaterials()
{
  /// Read gas properties from mother volume
  gas_         = mother_logic_->GetMaterial();
  pressure_    = gas_->GetPressure();
  temperature_ = gas_->GetTemperature();

  /// High density polyethylene for the field cage
  hdpe_ = materials::HDPE();

  /// PE500 for the holders
  pe500_ = materials::PE500();

  /// Copper for field rings
  copper_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");

  /// Teflon for the light tube
  teflon_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
  // teflon is the material used in the light-tube, and is covered by a G4LogicalSkinSurface
  // In Geant4 11.0.0, a bug in treating the OpBoundaryProcess produced in the surface makes the code fail.
  // This is avoided by setting an empty G4MaterialPropertiesTable of the G4Material.
  teflon_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());

  /// TPB coating
  tpb_ = materials::TPB();
  tpb_->SetMaterialPropertiesTable(opticalprops::TPB());

  /// Steel
  steel_ = materials::Steel316Ti();
  // In Geant4 11.0.0, a bug in treating the OpBoundaryProcess produced in the surface makes the code fail.
  // This is avoided by setting an empty G4MaterialPropertiesTable of the G4Material.
  steel_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
}


void Next100FieldCage::BuildActive()
{
  G4double new_active_zpos_ = z_fend - (active_length_/2. + fiber_end_z/2.);
  // G4double new_active_zpos_ = active_zpos_;

  /// Position of z planes
  G4double zplane[2] = {-active_length_/2. + gate_teflon_dist_ - overlap_,
                         active_length_/2.-(cathode_thickn_-grid_thickn_)/2. - fiber_end_z};
  // G4double zplane[2] = {-teflon_drift_length_/2.,
  //                        teflon_drift_length_/2.};
  /// Inner radius
  G4double rinner[2] = {0., 0.};
  /// Outer radius
  // G4double router[2] = {active_diam_/2., active_diam_/2.};
  G4double router[2] = {hh - fiber_diameter_/2.,
                        hh - fiber_diameter_/2.};
  // G4double router[2] = {active_diam_/2. - (teflon_thickn_ + fiber_diameter_),
  //                       active_diam_/2. - (teflon_thickn_ + fiber_diameter_)};

  G4Polyhedra* active_solid =
    new G4Polyhedra("ACTIVE_POLY", 0., twopi, n_panels, 2, zplane, rinner, router);

  G4Tubs* active_cathode_solid =
  new G4Tubs("ACT_CATHODE_RING", 0, cathode_int_diam_/2.,
              ((cathode_thickn_ - grid_thickn_)/2. + overlap_)/2., 0, twopi);

  G4ThreeVector act_cathode_pos =
  G4ThreeVector(0., 0., active_length_/2.-((cathode_thickn_ - grid_thickn_)/2.)/2. - overlap_/2.);

  G4UnionSolid* union_active =
    new G4UnionSolid ("ACTIVE", active_solid, active_cathode_solid, 0, act_cathode_pos);


// THIS ADENDA TO THE VOLUME IS OVERLAPPING WITH THE EL-GAP*************************************************************

  // //This volume is added as an extension of the active volume that reaches the gate grid.
  // G4Tubs* active_gate_solid =
  //   new G4Tubs("ACT_GATE_GAS", 0, gate_int_diam_/2., gate_teflon_dist_/2., 0, twopi);
  //
  // G4ThreeVector act_gate_pos =
  // G4ThreeVector(0., 0., -active_length_/2.+ gate_teflon_dist_/2.);
  //
  // union_active =
  //   new G4UnionSolid ("ACTIVE", union_active, active_gate_solid, 0, act_gate_pos);

// THIS ADENDA TO THE VOLUME IS OVERLAPPING WITH THE EL-GAP*************************************************************

  G4LogicalVolume* active_logic =
    new G4LogicalVolume(union_active, gas_, "ACTIVE");

  active_phys_ =
    new G4PVPlacement(0, G4ThreeVector(0., 0., new_active_zpos_),
                      active_logic, "ACTIVE", mother_logic_,
                      false, 0, false);

  /// Limit the step size in this volume for better tracking precision
  active_logic->SetUserLimits(new G4UserLimits(max_step_size_));


  /// Set the volume as an ionization sensitive detector
  IonizationSD* ionisd = new IonizationSD("/NEXT100/ACTIVE");
  active_logic->SetSensitiveDetector(ionisd);
  G4SDManager::GetSDMpointer()->AddNewDetector(ionisd);

  /// Define a drift field for this volume
  UniformElectricDriftField* field = new UniformElectricDriftField();
  G4double global_active_zpos = new_active_zpos_ - GetELzCoord();
  field->SetCathodePosition(global_active_zpos + active_length_/2.);
  field->SetAnodePosition(global_active_zpos - active_length_/2.);
  field->SetDriftVelocity(1. * mm/microsecond);
  field->SetTransverseDiffusion(drift_transv_diff_);
  field->SetLongitudinalDiffusion(drift_long_diff_);
  G4Region* drift_region = new G4Region("DRIFT");
  drift_region->SetUserInformation(field);
  drift_region->AddRootLogicalVolume(active_logic);


  /// Vertex generator
  active_gen_ = new CylinderPointSampler2020(0., active_diam_/2., active_length_/2.,
                                             0., twopi, nullptr,
                                             G4ThreeVector(0., 0., new_active_zpos_));


  /// Visibilities
  // active_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  active_logic->SetVisAttributes(nexus::Yellow());

  /// Verbosity
  if (verbosity_) {
    G4cout << "Active starts in " << (new_active_zpos_ - active_length_/2.)/mm
           << " mm and ends in "
           << (new_active_zpos_ + active_length_/2.)/mm << " mm" << G4endl;
  }
}


void Next100FieldCage::BuildCathode()
{
  G4Tubs* cathode_solid =
    new G4Tubs("CATHODE_RING", cathode_int_diam_/2.,cathode_ext_diam_/2.,
               cathode_thickn_/2., 0, twopi);

  G4LogicalVolume* cathode_logic =
    new G4LogicalVolume(cathode_solid, steel_, "CATHODE_RING");

  new G4PVPlacement(0, G4ThreeVector(0., 0., cathode_zpos_),
                    cathode_logic, "CATHODE_RING", mother_logic_,
                    false, 0, false);

  G4Material* fgrid_mat = materials::FakeDielectric(gas_, "cath_grid_mat");
  fgrid_mat->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_, temperature_,
                                                               cath_grid_transparency_,
                                                               grid_thickn_));

  G4Tubs* diel_grid_solid =
    new G4Tubs("CATHODE_GRID", 0., cathode_int_diam_/2., grid_thickn_/2., 0, twopi);

  G4LogicalVolume* diel_grid_logic =
    new G4LogicalVolume(diel_grid_solid, fgrid_mat, "CATHODE_GRID");

  new G4PVPlacement(0, G4ThreeVector(0., 0., cathode_zpos_),
                    diel_grid_logic, "CATHODE_GRID", mother_logic_,
                    false, 0, false);

  // Cathode ring vertex generator
  cathode_gen_ = new CylinderPointSampler2020(cathode_int_diam_/2.,cathode_ext_diam_/2.,
                                           cathode_thickn_/2.,0., twopi, nullptr,
                                           G4ThreeVector(0., 0., cathode_zpos_));


  /// Visibilities
  if (visibility_) {
    G4VisAttributes grey = nexus::LightGrey();
    G4VisAttributes cathode_col = nexus::DarkGrey();
    cathode_col.SetForceSolid(true);
    diel_grid_logic->SetVisAttributes(grey);
    cathode_logic->SetVisAttributes(cathode_col);
  } else {
    G4VisAttributes cathode_col = nexus::DarkGrey();
    cathode_col.SetForceSolid(true);
    diel_grid_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    cathode_logic->SetVisAttributes(cathode_col);
  }


  /// Verbosity
  if (verbosity_) {
    G4cout << "Cathode grid pos z: " << (cathode_zpos_)/mm << " mm" << G4endl;
  }

}


void Next100FieldCage::BuildBuffer()
{
  G4double buffer_zpos = active_zpos_ + active_length_/2. + grid_thickn_ + buffer_length_/2.;

  /// Position of z planes
  G4double zplane[2] = {-buffer_length_/2.+(cathode_thickn_-grid_thickn_)/2., buffer_length_/2.};
  /// Inner radius
  G4double rinner[2] = {0., 0.};
  /// Outer radius
  G4double router[2] = {active_diam_/2., active_diam_/2.};

  G4Polyhedra* buffer_solid =
    new G4Polyhedra("BUFFER_POLY", 0., twopi, n_panels, 2, zplane, rinner, router);

  G4Tubs* buffer_cathode_solid =
    new G4Tubs("BUFF_CATHODE_RING", 0, cathode_int_diam_/2.,
              (cathode_thickn_/2. - grid_thickn_/2.)/2. +  overlap_/2., 0, twopi);

  G4ThreeVector buff_cathode_pos =
  G4ThreeVector(0., 0., -buffer_length_/2. + (cathode_thickn_/2.-grid_thickn_/2.)/2. +overlap_/2.);

  G4UnionSolid* union_buffer =
    new G4UnionSolid("BUFFER", buffer_solid, buffer_cathode_solid, 0, buff_cathode_pos);

  G4LogicalVolume* buffer_logic =
    new G4LogicalVolume(union_buffer, gas_, "BUFFER");

  buffer_phys_ =
    new G4PVPlacement(0, G4ThreeVector(0., 0., buffer_zpos),
                      buffer_logic, "BUFFER", mother_logic_,
                      false, 0, false);

  /// Set the volume as an ionization sensitive detector
  IonizationSD* buffsd = new IonizationSD("/NEXT100/BUFFER");
  buffsd->IncludeInTotalEnergyDeposit(false);
  buffer_logic->SetSensitiveDetector(buffsd);
  G4SDManager::GetSDMpointer()->AddNewDetector(buffsd);


  /// Vertex generator
  G4double active_ext_radius = active_diam_/2. / cos(pi/n_panels);
  buffer_gen_ = new CylinderPointSampler2020(0., active_ext_radius, buffer_length_/2.,
                                             0., twopi, nullptr,
                                             G4ThreeVector(0., 0., buffer_zpos));

  /// Vertex generator for all xenon
  G4double xenon_length = el_gap_length_ + active_length_ +
                          grid_thickn_ + buffer_length_ ;
  G4double xenon_zpos   = (el_gap_length_ * el_gap_zpos_ +
                          active_length_ * active_zpos_ +
                          grid_thickn_ * cathode_zpos_ +
                          buffer_length_ * buffer_zpos) / xenon_length;
  xenon_gen_ = new CylinderPointSampler2020(0., active_ext_radius, xenon_length,
                                            0., twopi, nullptr,
                                            G4ThreeVector(0., 0., xenon_zpos));

  /// Visibilities
  buffer_logic->SetVisAttributes(G4VisAttributes::GetInvisible());

  /// Verbosity
  if (verbosity_) {
    G4cout << "Buffer (gas) starts in " << buffer_zpos - buffer_length_/2.
           << " and ends in "
           << buffer_zpos + buffer_length_/2. << G4endl;
  }
}


void Next100FieldCage::BuildELRegion()
{
  /// GATE ring.
  G4Tubs* gate_solid =
    new G4Tubs("GATE_RING", gate_int_diam_/2., gate_ext_diam_/2., gate_ring_thickn_/2., 0, twopi);

  G4LogicalVolume* gate_logic =
    new G4LogicalVolume(gate_solid, steel_, "GATE_RING");

  new G4PVPlacement(0, G4ThreeVector(0., 0., gate_zpos_),
                    gate_logic, "GATE_RING", mother_logic_,
                    false, 0, false);

  /// EL gap.
  G4Tubs* el_gap_solid =
    new G4Tubs("EL_GAP", 0., gate_int_diam_/2., (el_gap_length_ + 2*grid_thickn_)/2., 0, twopi);

  G4LogicalVolume* el_gap_logic =
    new G4LogicalVolume(el_gap_solid, gas_, "EL_GAP");

  new G4PVPlacement(0, G4ThreeVector(0., 0., el_gap_zpos_),
                    el_gap_logic, "EL_GAP", mother_logic_,
                    false, 0, false);

  /// ANODE ring.
  G4Tubs* anode_solid =
    new G4Tubs("ANODE_RING", gate_int_diam_/2., gate_ext_diam_/2., gate_ring_thickn_/2., 0, twopi);

  G4LogicalVolume* anode_logic =
    new G4LogicalVolume(anode_solid, steel_, "ANODE_RING");

  new G4PVPlacement(0, G4ThreeVector(0., 0., anode_zpos_),
                    anode_logic, "ANODE_RING", mother_logic_,
                    false, 0, false);

  if (elfield_) {
    /// ma EL electric field
    UniformElectricDriftField* el_field = new UniformElectricDriftField();
    G4double global_el_gap_zpos = el_gap_zpos_ - GetELzCoord();
    el_field->SetCathodePosition(global_el_gap_zpos + el_gap_length_/2. + grid_thickn_);
    el_field->SetAnodePosition  (global_el_gap_zpos - el_gap_length_/2. - grid_thickn_);
    el_field->SetDriftVelocity(2.5 * mm/microsecond);
    el_field->SetTransverseDiffusion(ELtransv_diff_);
    el_field->SetLongitudinalDiffusion(ELlong_diff_);
    el_field->SetLightYield(XenonELLightYield(ELelectric_field_, pressure_));
    G4Region* el_region = new G4Region("EL_REGION");
    el_region->SetUserInformation(el_field);
    el_region->AddRootLogicalVolume(el_gap_logic);
  }

  /// EL grids
  G4Material* fgrid_mat = materials::FakeDielectric(gas_, "el_grid_mat");
  fgrid_mat->SetMaterialPropertiesTable(opticalprops::FakeGrid(pressure_, temperature_,
                                                               el_grid_transparency_,
                                                               grid_thickn_));

  /// Dimensions & position: the grids are simulated inside the EL gap.
  /// Their thickness is symbolic.
  G4Tubs* diel_grid_solid =
    new G4Tubs("EL_GRID", 0., gate_int_diam_/2., grid_thickn_/2., 0, twopi);

  G4LogicalVolume* diel_grid_logic =
    new G4LogicalVolume(diel_grid_solid, fgrid_mat, "EL_GRID");

  new G4PVPlacement(0, G4ThreeVector(0., 0., el_gap_length_/2. + grid_thickn_/2.),
                    diel_grid_logic, "EL_GRID_GATE", el_gap_logic,
                    false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(0., 0., -el_gap_length_/2. - grid_thickn_/2.),
                    diel_grid_logic, "EL_GRID_ANODE", el_gap_logic,
                    false, 1, false);

  // Vertex generator
  if (el_gap_gen_disk_zmin_ > el_gap_gen_disk_zmax_)
    G4Exception("[Next100FieldCage]", "Next100FieldCage()", FatalErrorInArgument,
                "Error in configuration of EL gap generator: zmax < zmin");

  G4double el_gap_gen_disk_thickn = el_gap_length_ *
                                    (el_gap_gen_disk_zmax_ - el_gap_gen_disk_zmin_);

  G4double el_gap_gen_disk_z = el_gap_zpos_ + el_gap_length_/2.-
                               el_gap_length_ * el_gap_gen_disk_zmin_ -
                               el_gap_gen_disk_thickn/2.;

  G4ThreeVector el_gap_gen_pos(el_gap_gen_disk_x_, el_gap_gen_disk_y_, el_gap_gen_disk_z);

  el_gap_gen_ = new CylinderPointSampler2020(0., el_gap_gen_disk_diam_/2.,
                                             el_gap_gen_disk_thickn/2., 0., twopi,
                                             nullptr, el_gap_gen_pos);

  // Gate ring vertex generator
  gate_gen_ = new CylinderPointSampler2020(gate_int_diam_/2., gate_ext_diam_/2., gate_ring_thickn_/2.,
                                           0., twopi, nullptr, G4ThreeVector(0., 0., gate_zpos_));
  // Anode ring vertex generator
  anode_gen_ = new CylinderPointSampler2020(gate_int_diam_/2., gate_ext_diam_/2., gate_ring_thickn_/2.,
                                            0., twopi, nullptr, G4ThreeVector(0., 0., anode_zpos_));

  /// Visibilities
  if (visibility_) {
    G4VisAttributes light_blue = nexus::LightBlue();
    el_gap_logic->SetVisAttributes(light_blue);
    diel_grid_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  } else {
    el_gap_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    diel_grid_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  }
  G4VisAttributes grey = nexus::DarkGrey();
  grey.SetForceSolid(true);
  gate_logic->SetVisAttributes(grey);
  anode_logic->SetVisAttributes(grey);

  /// Verbosity
  if (verbosity_) {
    G4cout << "EL gap starts in " << (el_gap_zpos_ - el_gap_length_/2.)/mm
           << " mm and ends in " << (el_gap_zpos_ + el_gap_length_/2.)/mm << G4endl;
    G4cout << G4endl;
  }
}


void Next100FieldCage::BuildFiberBarrel()
{

    // MATERIALS /////////////////////////////////////////////
    G4OpticalSurface* opsur_teflon =
      new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
    opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());

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

   // TEFLON PANEL /////////////////////////////////////////////
   G4Box* teflon_panel =
     new G4Box("TEFLON_PANEL", panel_width_/2., panel_length_/2., panel_thickness_/2.);
   G4LogicalVolume* teflon_panel_logic =
     new G4LogicalVolume(teflon_panel, teflon_, "TEFLON_PANEL");

   new G4LogicalSkinSurface("TEFLON_PANEL_OPSURF", teflon_panel_logic, opsur_teflon);

   teflon_panel_logic->SetVisAttributes(nexus::White());
   // teflon_panel_logic->SetVisAttributes(G4VisAttributes::GetInvisible());


    // DETECTOR ////////////////////////////////////////////////////////////////
    G4double sensor_width = panel_width_/n_sensors;
    G4String sensor_name = "F_SENSOR";

    /// Build the sensor
    photo_sensor_  = new GenericPhotosensor(sensor_name, sensor_width, fiber_diameter_, sens_z);


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
      G4double t_binning = .1 * ns;
      // G4double t_binning = 100. * ns;
      // G4double t_binning = 1. * ns;

      photo_sensor_ ->SetTimeBinning(t_binning); // Size of fiber sensors time binning

      G4cout << "[FiberBarrel] Using " << t_binning << " [ns] binning"<< G4endl;


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

    G4Tubs* fiber_end_solid_vol =
      new G4Tubs("fiber_end", 0, fiber_diameter_ / 2, fiber_end_z/2., 0, 2 * M_PI);

    G4LogicalVolume* fiber_end_logic_vol =
      new G4LogicalVolume(fiber_end_solid_vol, fiber_end_mat, "FIBER_END");
    G4OpticalSurface* opsur_al =
      new G4OpticalSurface("POLISHED_AL_OPSURF", unified, polished, dielectric_metal);
    opsur_al->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("POLISHED_AL_OPSURF", fiber_end_logic_vol, opsur_al);

    fiber_end_logic_vol  ->SetVisAttributes(nexus::Blue());


    // fiber ////////////////////////////////////////////////////

    fiber_ = new GenericWLSFiber(fiber_type_, true, fiber_diameter_, fiber_length, true, coated_, this_coating, this_fiber, true);

    fiber_->SetCoreOpticalProperties(this_fiber_optical);
    fiber_->SetCoatingOpticalProperties(this_coating_optical);

    fiber_->Construct();
    G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
    if (fiber_type_ == "Y11")
      fiber_logic->SetVisAttributes(nexus::LightGreenAlpha());
    else if (fiber_type_ == "B2")
      fiber_logic->SetVisAttributes(nexus::LightBlueAlpha());


   // teflon cap to cover EP ////////////////////////////////////////////
   G4Tubs* teflon_cap =
   new G4Tubs("TEFLON_CAP", 0, hh - fiber_diameter_/2., fiber_end_z/2., 0, twopi);
   G4Material* teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
   teflon->SetMaterialPropertiesTable(opticalprops::PTFE());
   G4LogicalVolume* teflon_cap_logic =
     new G4LogicalVolume(teflon_cap, teflon_, "TEFLON");

  //  G4OpticalSurface* opsur_teflon =
  //    new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
  //  opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());

   new G4LogicalSkinSurface("TEFLON_OPSURF", teflon_cap_logic, opsur_teflon);

   if (cap_visibility_ == false)
   {
     teflon_cap_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
   }

   new G4PVPlacement(0, G4ThreeVector(0, 0, z_fend),
                     teflon_cap_logic, "TEFLON_CAP", mother_logic_,
                     true, 0, true);


   // PLACEMENT /////////////////////////////////////////////
   G4double rot_angle;
   G4double theta0 =  (10.)*pi/180.;

   for (G4int itheta=0; itheta < n_panels; itheta++) {
   // for (G4int itheta=0; itheta < 3; itheta++) {

     // panels
     G4double theta = theta0 + dif_theta * itheta;
     G4double x = h * std::cos(theta) * mm;
     G4double y = h * std::sin(theta) * mm;
     G4double phi = pi/2. + std::atan2(y, x);
     std::string label = std::to_string(itheta);

     G4RotationMatrix* panel_rot = new G4RotationMatrix();
     rot_angle = pi/2.;
     panel_rot->rotateX(rot_angle);
     panel_rot->rotateY(phi);
     new G4PVPlacement(panel_rot, G4ThreeVector(x,y, z),
                       teflon_panel_logic, "PANEL-"+label, mother_logic_,
                       false, itheta, false);

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

         new G4PVPlacement(0, G4ThreeVector(xx_f, yy_f, z_f),
                           fiber_logic, "FIBER-" + label + label2, mother_logic_,
                           false, n_panels + ii, false);
        new G4PVPlacement(0, G4ThreeVector(xx_f, yy_f, z_fend),
                          fiber_end_logic_vol, "ALUMINIUMR-" + label + label2, mother_logic_,
                          false, n_panels + n_fibers*itheta + ii, false);
       }
     for (G4int jj=0; jj < n_sensors; jj++) {
    // for (G4int jj=0; jj < 3; jj++) {

          G4double xx_s = x0_s - dl_sens*jj*std::cos(phi);
          G4double yy_s = y0_s - dl_sens*jj*std::sin(phi);

          std::string label3 = std::to_string(jj);

          G4RotationMatrix* sensor_rot = new G4RotationMatrix();
          // rot_angle = 0.;
          rot_angle = M_PI;
          sensor_rot->rotateY(rot_angle);
          sensor_rot->rotateZ(phi);
          new G4PVPlacement(sensor_rot, G4ThreeVector(xx_s, yy_s, z_s),
                            photo_sensor_logic, "SENS-" + label + label3, mother_logic_,
                            true, n_panels*(1 + n_fibers) + n_sensors*itheta  + jj, true);

    }

   }

}
void Next100FieldCage::BuildLightTube()
{
  /// DRIFT PART ///
  /// Position of z planes
  G4double zplane[2] = {-teflon_drift_length_/2., teflon_drift_length_/2.};
  /// Inner radius
  G4double rinner[2] = {active_diam_/2., active_diam_/2.};
  /// Outer radius
  G4double router[2] =
    {(active_diam_ + 2.*teflon_thickn_)/2., (active_diam_ + 2.*teflon_thickn_)/2.};

  G4Polyhedra* teflon_drift_solid =
    new G4Polyhedra("LIGHT_TUBE_DRIFT", 0., twopi, n_panels, 2, zplane, rinner, router);

  G4LogicalVolume* teflon_drift_logic =
    new G4LogicalVolume(teflon_drift_solid, teflon_, "LIGHT_TUBE_DRIFT");

  new G4PVPlacement(0, G4ThreeVector(0., 0., teflon_drift_zpos_),
                    teflon_drift_logic, "LIGHT_TUBE_DRIFT", mother_logic_,
                    false, 0, false);


  /// TPB on teflon surface
  G4double router_tpb[2] =
    {(active_diam_ + 2.*tpb_thickn_)/2., (active_diam_ + 2.*tpb_thickn_)/2.};

  G4Polyhedra* tpb_drift_solid =
    new  G4Polyhedra("DRIFT_TPB", 0., twopi, n_panels, 2, zplane, rinner, router_tpb);
  G4LogicalVolume* tpb_drift_logic =
    new G4LogicalVolume(tpb_drift_solid, tpb_, "DRIFT_TPB");
  G4VPhysicalVolume* tpb_drift_phys =
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                      tpb_drift_logic, "DRIFT_TPB", teflon_drift_logic,
                      false, 0, false);

  /// BUFFER PART ///
  G4double zplane_buff[2] = {-teflon_buffer_length_/2., teflon_buffer_length_/2.};
  G4double router_buff[2] =
    {(active_diam_ + 2.*teflon_thickn_)/2., (active_diam_ + 2.*teflon_thickn_)/2.};

  G4Polyhedra* teflon_buffer_solid =
   new G4Polyhedra("LIGHT_TUBE_BUFFER", 0., twopi, n_panels, 2, zplane_buff, rinner, router_buff);

  G4LogicalVolume* teflon_buffer_logic =
    new G4LogicalVolume(teflon_buffer_solid, teflon_, "LIGHT_TUBE_BUFFER");

  new G4PVPlacement(0, G4ThreeVector(0., 0., teflon_buffer_zpos_),
                    teflon_buffer_logic, "LIGHT_TUBE_BUFFER", mother_logic_,
                    false, 0, false);

  /// TPB on teflon surface
  G4double router_tpb_buff[2] =
    {(active_diam_ + 2.*tpb_thickn_)/2., (active_diam_ + 2.*tpb_thickn_)/2.};

  G4Polyhedra* tpb_buffer_solid =
    new  G4Polyhedra("BUFFER_TPB", 0., twopi, n_panels, 2, zplane_buff, rinner, router_tpb_buff);

  G4LogicalVolume* tpb_buffer_logic =
    new G4LogicalVolume(tpb_buffer_solid, tpb_, "BUFFER_TPB");

  G4VPhysicalVolume* tpb_buffer_phys =
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                      tpb_buffer_logic, "BUFFER_TPB", teflon_buffer_logic,
                      false, 0, false);

  /// Optical surface on teflon ///
  G4OpticalSurface* refl_Surf =
    new G4OpticalSurface("refl_Surf", unified, ground, dielectric_metal, .01);
  refl_Surf->SetMaterialPropertiesTable(opticalprops::PTFE());
  new G4LogicalSkinSurface("refl_teflon_surf", teflon_drift_logic, refl_Surf);
  new G4LogicalSkinSurface("refl_teflon_surf", teflon_buffer_logic, refl_Surf);

  /// Optical surface between xenon and TPB to model roughness ///
  G4OpticalSurface* gas_tpb_teflon_surf =
    new G4OpticalSurface("gas_tpb_teflon_surf", glisur, ground,
                         dielectric_dielectric, .01);

  new G4LogicalBorderSurface("gas_tpb_teflon_surf", tpb_drift_phys, active_phys_,
                             gas_tpb_teflon_surf);
  new G4LogicalBorderSurface("gas_tpb_teflon_surf", active_phys_, tpb_drift_phys,
                             gas_tpb_teflon_surf);
  new G4LogicalBorderSurface("gas_tpb_teflon_surf", tpb_buffer_phys, buffer_phys_,
                             gas_tpb_teflon_surf);
  new G4LogicalBorderSurface("gas_tpb_teflon_surf", buffer_phys_, tpb_buffer_phys,
                             gas_tpb_teflon_surf);

  // Vertex generator
  G4double teflon_ext_radius = (active_diam_ + 2.*teflon_thickn_)/2. / cos(pi/n_panels);
  G4double cathode_gap_zpos  = teflon_drift_zpos_ + teflon_drift_length_/2. + cathode_thickn_/2.;
  G4double teflon_zpos = (teflon_drift_length_ * teflon_drift_zpos_ +
                         cathode_thickn_ * cathode_gap_zpos +
                         teflon_buffer_length_ * teflon_buffer_zpos_) / teflon_total_length_;

  teflon_gen_ = new CylinderPointSampler2020(active_diam_/2., teflon_ext_radius,
                                             teflon_total_length_/2., 0., twopi,
                                             nullptr, G4ThreeVector (0., 0., teflon_zpos));

  // Visibilities
  if (visibility_) {
    G4VisAttributes light_yellow = nexus::YellowAlpha();
    light_yellow.SetForceSolid(true);
    teflon_drift_logic->SetVisAttributes(light_yellow);
    teflon_buffer_logic->SetVisAttributes(light_yellow);
    tpb_drift_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    tpb_buffer_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  }
  else {
    teflon_drift_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    teflon_buffer_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    tpb_drift_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    tpb_buffer_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  }

}


void Next100FieldCage::BuildFieldCage()
{
  // HDPE cylinder.
  G4double hdpe_tube_z_pos = teflon_buffer_zpos_ - (hdpe_length_ - teflon_buffer_length_)/2.;

  G4Tubs* hdpe_tube_solid =
    new G4Tubs("HDPE_TUBE", hdpe_tube_int_diam_/2., hdpe_tube_ext_diam_/2.,
               hdpe_length_/2., 0, twopi);

  G4LogicalVolume* hdpe_tube_logic =
    new G4LogicalVolume(hdpe_tube_solid, hdpe_, "HDPE_TUBE");

  new G4PVPlacement(0, G4ThreeVector(0., 0., hdpe_tube_z_pos),
                    hdpe_tube_logic, "HDPE_TUBE", mother_logic_,
                    false, 0, false);

  hdpe_gen_ = new CylinderPointSampler2020(hdpe_tube_int_diam_/2., hdpe_tube_ext_diam_/2.,
                                           hdpe_length_/2.,0., twopi, nullptr,
                                           G4ThreeVector(0., 0., hdpe_tube_z_pos));

  G4double active_short_z = 13.5 * mm; //Thickness of holder first holder in the active volume.
  G4double buffer_short_z = 37.  * mm;
  G4double ring_drift_buffer_dist = 72.*mm;
  G4int    num_drift_rings = 48;
  G4int    num_buffer_rings = 4;
  G4double posz;
  G4double first_ring_drift_z_pos = GetELzCoord() + gate_teflon_dist_ + drift_ring_dist_/2. + active_short_z/2.;

  G4double first_ring_buff_z_pos = first_ring_drift_z_pos + (num_drift_rings-1)*drift_ring_dist_ +
                                   ring_drift_buffer_dist;

  G4Tubs* ring_solid =
    new G4Tubs("FIELD_RING", ring_int_diam_/2., ring_ext_diam_/2., ring_thickn_/2., 0, twopi);

  G4LogicalVolume* ring_logic =
    new G4LogicalVolume(ring_solid, copper_, "FIELD_RING");

  //Placement of the drift rings.
  for (G4int i=0; i<num_drift_rings; i++) {
    posz = first_ring_drift_z_pos + i*drift_ring_dist_;
    new G4PVPlacement(0, G4ThreeVector(0., 0., posz),
                      ring_logic, "FIELD_RING", mother_logic_,
                      false, i, false);
  }

  //Placement of the buffer rings.
  for (G4int i=0; i<num_buffer_rings; i++) {
    posz = first_ring_buff_z_pos + i*buffer_ring_dist_;
    new G4PVPlacement(0, G4ThreeVector(0., 0., posz),
                      ring_logic, "FIELD_RING", mother_logic_,
                      false, i, false);
  }

  // ring vertex generator
  G4double ring_gen_lenght =   first_ring_buff_z_pos + (num_buffer_rings-1)*buffer_ring_dist_
                             - first_ring_drift_z_pos + ring_thickn_;
  G4double ring_gen_zpos = first_ring_drift_z_pos + ring_gen_lenght/2. - ring_thickn_/2.;
  ring_gen_ = new CylinderPointSampler2020(ring_int_diam_/2., ring_ext_diam_/2., ring_gen_lenght/2.,
                                           0., twopi, nullptr,
                                           G4ThreeVector(0., 0., ring_gen_zpos));

  // Ring holders.
  // ACTIVE holders.
  G4Box* active_short_solid =
    new G4Box("ACT_SHORT", holder_x_/2., holder_short_y_/2.+overlap_/2., active_short_z/2.);

  G4double first_act_short_z = -teflon_drift_length_/2.+ active_short_z/2.;

  G4Box* active_long_solid =
    new G4Box("ACT_LONG", holder_x_/2., holder_long_y_/2., teflon_drift_length_/2.);

  G4UnionSolid* act_holder_solid =
    new G4UnionSolid ("ACT_HOLDER", active_long_solid, active_short_solid, 0,
                      G4ThreeVector(0.,holder_long_y_/2.+holder_short_y_/2.-overlap_/2.,
                                    first_act_short_z));

  for (G4int j=1; j<num_drift_rings; j++) {
    posz = first_act_short_z + j*drift_ring_dist_;

    act_holder_solid =
      new G4UnionSolid("ACT_HOLDER", act_holder_solid, active_short_solid, 0,
                       G4ThreeVector(0.,holder_long_y_/2.+holder_short_y_/2.-overlap_/2., posz));
    }

  G4LogicalVolume* act_holder_logic =
    new G4LogicalVolume(act_holder_solid, pe500_, "ACT_HOLDER");
  G4int numbering=0;
  for (G4int i=10; i<360; i +=20){
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot -> rotateZ((90-i) *deg);
    new G4PVPlacement(rot,G4ThreeVector(holder_r_*cos(i*deg),holder_r_*sin(i*deg),teflon_drift_zpos_),
                      act_holder_logic, "ACT_HOLDER", mother_logic_,
                      false, numbering, false);
    numbering +=1;}

  // BUFFER holders.
  G4Box* buffer_short_solid =
    new G4Box("BUFF_SHORT", holder_x_/2., holder_short_y_/2.+overlap_/2., buffer_short_z/2.);

  G4double first_buff_short_z = -teflon_buffer_length_/2. +
                                (ring_drift_buffer_dist/2.-cathode_thickn_/2.) +
                                buffer_ring_dist_/2.;
  G4Box* buffer_long_solid =
    new G4Box("BUFF_LONG", holder_x_/2., holder_long_y_/2., teflon_buffer_length_/2.);

  G4UnionSolid* buff_holder_solid =
    new G4UnionSolid ("BUFF_HOLDER", buffer_long_solid, buffer_short_solid, 0,
                      G4ThreeVector(0.,holder_long_y_/2.+holder_short_y_/2.-overlap_/2.,
                                    first_buff_short_z));

  for (G4int j=1; j<num_buffer_rings-1; j++) {
    posz = first_buff_short_z + j*buffer_ring_dist_;

    buff_holder_solid =
      new G4UnionSolid("BUFF_HOLDER", buff_holder_solid, buffer_short_solid, 0,
                       G4ThreeVector(0.,holder_long_y_/2.+holder_short_y_/2.-overlap_/2.,posz));
    }

  G4double buffer_last_z  = 63.2 *mm;
  G4Box* buffer_last_solid =
    new G4Box("BUFF_LAST", holder_x_/2., holder_short_y_/2.+overlap_/2., buffer_last_z/2.);

  buff_holder_solid =
    new G4UnionSolid("BUFF_HOLDER", buff_holder_solid, buffer_last_solid, 0,
                     G4ThreeVector(0.,holder_long_y_/2. + holder_short_y_/2.-overlap_/2.,
                                   teflon_buffer_length_/2. - buffer_last_z/2.));

  G4LogicalVolume* buff_holder_logic =
    new G4LogicalVolume(buff_holder_solid, pe500_, "BUFF_HOLDER");

  numbering=0;
  for (G4int i=10; i<360; i +=20){
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot -> rotateZ((90-i) *deg);
    new G4PVPlacement(rot, G4ThreeVector(holder_r_*cos(i*deg), holder_r_*sin(i*deg),
                      teflon_buffer_zpos_),buff_holder_logic, "BUFF_HOLDER", mother_logic_,
                      false, numbering, false);
    numbering +=1;}

  // CATHODE holders.
  G4double cathode_long_y = 29.*mm;
  G4double cathode_long_z = 61*mm;
  G4double cathode_short_z = 24.5*mm;
  G4Box* cathode_large_solid =
    new G4Box("CATHODE_LARGE", holder_x_/2., cathode_long_y/2., cathode_long_z/2.);

  G4Box* cathode_short_solid =
    new G4Box("CATHODE_SHORT", holder_x_/2., holder_short_y_/2., cathode_short_z/2.);

  G4UnionSolid* cathode_holder_solid =
    new G4UnionSolid ("CATHODE_HOLDER", cathode_large_solid, cathode_short_solid, 0,
                      G4ThreeVector(0.,-(holder_short_y_/2.-cathode_long_y/2),
                                    cathode_long_z/2.-cathode_short_z/2.));

  cathode_holder_solid =
    new G4UnionSolid("CATHODE_HOLDER", cathode_holder_solid, cathode_short_solid, 0,
                      G4ThreeVector(0.,-(holder_short_y_/2.-cathode_long_y/2),
                                    -(cathode_long_z/2.-cathode_short_z/2.)));

  G4LogicalVolume* cathode_holder_logic =
    new G4LogicalVolume(cathode_holder_solid, pe500_, "CATHODE_HOLDER");

  numbering=0;
  G4double cathode_holder_r = (active_diam_+2*teflon_thickn_+ 2*holder_long_y_+
                              2*holder_short_y_)/2.-cathode_long_y/2.;
  for (G4int i=10; i<360; i +=20){
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot -> rotateZ((90-i) *deg);
    new G4PVPlacement(rot, G4ThreeVector(cathode_holder_r*cos(i*deg),cathode_holder_r*sin(i*deg),
                      cathode_zpos_),cathode_holder_logic, "CATHODE_HOLDER", mother_logic_,
                      false, numbering, false);
    numbering +=1;}

  holder_gen_ = new CylinderPointSampler2020(holder_r_ - holder_long_y_/2.,
                                             holder_r_ + holder_long_y_/2. + holder_short_y_,
                                             gate_sapphire_wdw_dist_/2., 0., twopi, nullptr,
                                             G4ThreeVector(0., 0., gate_grid_zpos_ + gate_sapphire_wdw_dist_/2.));

  /// Visibilities
  if (visibility_) {
    G4VisAttributes ring_col = nexus::CopperBrown();
    ring_col.SetForceSolid(true);
    ring_logic->SetVisAttributes(ring_col);
    G4VisAttributes hdpe_col =nexus::WhiteAlpha();
    hdpe_col.SetForceSolid(true);
    hdpe_tube_logic->SetVisAttributes(hdpe_col);
    G4VisAttributes hold_col = nexus::LightGrey();
    hold_col.SetForceSolid(true);
    act_holder_logic->SetVisAttributes(hold_col);
    buff_holder_logic->SetVisAttributes(hold_col);
    cathode_holder_logic->SetVisAttributes(hold_col);
  } else {
    ring_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    hdpe_tube_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    act_holder_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    buff_holder_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    cathode_holder_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
  }
}


Next100FieldCage::~Next100FieldCage()
{
  delete active_gen_;
  delete buffer_gen_;
  delete xenon_gen_;
  delete teflon_gen_;
  delete el_gap_gen_;
  delete hdpe_gen_;
  delete ring_gen_;
  delete cathode_gen_;
  delete gate_gen_;
  delete anode_gen_;
  delete holder_gen_;
}


G4ThreeVector Next100FieldCage::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex(0., 0., 0.);

  if (region == "CENTER") {
    vertex = G4ThreeVector(0., 0., active_zpos_);
  }

  else if (region == "ACTIVE") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = active_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "CATHODE_RING") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = cathode_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "BUFFER") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = buffer_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "XENON") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = xenon_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (
    VertexVolume->GetName() != "ACTIVE" &&
    VertexVolume->GetName() != "BUFFER" &&
    VertexVolume->GetName() != "EL_GAP");
  }

  else if (region == "LIGHT_TUBE") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = teflon_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (
    VertexVolume->GetName() != "LIGHT_TUBE_DRIFT" &&
    VertexVolume->GetName() != "LIGHT_TUBE_BUFFER" );
  }

  else if (region == "HDPE_TUBE") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = hdpe_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "EL_GAP") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = el_gap_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "FIELD_RING") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = ring_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "GATE_RING") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = gate_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "ANODE_RING") {
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = anode_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while (VertexVolume->GetName() != region);
  }

  else if (region == "RING_HOLDER"){
    G4VPhysicalVolume *VertexVolume;
    do {
      vertex = holder_gen_->GenerateVertex("VOLUME");
      G4ThreeVector glob_vtx(vertex);
      glob_vtx = glob_vtx + G4ThreeVector(0, 0, -GetELzCoord());
      VertexVolume =
        geom_navigator_->LocateGlobalPointAndSetup(glob_vtx, 0, false);
    } while ((VertexVolume->GetName() != "ACT_HOLDER")  &&
             (VertexVolume->GetName() != "BUFF_HOLDER") &&
             (VertexVolume->GetName() != "CATHODE_HOLDER"));
 }

  else {
    G4Exception("[Next100FieldCage]", "GenerateVertex()", FatalException,
    "Unknown vertex generation region!");
  }

  return vertex;
}


G4ThreeVector Next100FieldCage::GetActivePosition() const
{
  return G4ThreeVector (0., 0., active_zpos_);
}
