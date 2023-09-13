// ----------------------------------------------------------------------------
// nexus | Next100.cc
//
// Main class that constructs the geometry of the NEXT-100 detector.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100.h"
#include "BoxPointSampler.h"
#include "MuonsPointSampler.h"
#include "LSCHallA.h"
#include "Next100Shielding.h"
#include "Next100Vessel.h"
#include "Next100Ics.h"
#include "Next100InnerElements.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4UserLimits.hh>

// // Marian's adenda
// #include "OpticalMaterialProperties.h"
// #include "Visibilities.h"
//
// #include <G4OpticalSurface.hh>
// #include <G4LogicalSkinSurface.hh>
// //


namespace nexus {

  REGISTER_CLASS(Next100, GeometryBase)

  using namespace CLHEP;

  Next100::Next100():
    GeometryBase(),
    // Lab dimensions
    lab_size_ (5. * m),

    // common used variables in geomety components
    // 0.1 mm grid thickness
    // note that if grid thickness change it must be also changed in Next100FieldCage.cc
    gate_tracking_plane_distance_((26.1 + 0.1)   * mm),
    gate_sapphire_wdw_distance_  ((1458.2 - 0.1) * mm),

    specific_vertex_{},
    lab_walls_(false)//,

    // fiber_type_ ("Y11")
  {

    msg_ = new G4GenericMessenger(this, "/Geometry/Next100/",
				  "Control commands of geometry Next100.");

    msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_,
      "Set generation vertex.");

    msg_->DeclareProperty("lab_walls", lab_walls_, "Placement of Hall A walls");

    // msg_->DeclareProperty("fiber_type", fiber_type_, "Fiber type (Y11 or B2)");

  // The following methods must be invoked in this particular
  // order since some of them depend on the previous ones

  // Shielding
  shielding_ = new Next100Shielding();

  //Lab walls
  hallA_walls_ = new LSCHallA();

  // Vessel
  vessel_ = new Next100Vessel();

  // Internal copper shielding
  ics_ = new Next100Ics();

  // Inner Elements
  inner_elements_ = new Next100InnerElements();

  }


  Next100::~Next100()
  {
    delete inner_elements_;
    delete ics_;
    delete vessel_;
    delete shielding_;
    delete lab_gen_;
    delete hallA_walls_;
  }


  void Next100::Construct()
  {
     G4cout << "[Next100] *** Full Next100 simulation with fibers ***" << G4endl;
    // LAB /////////////////////////////////////////////////////////////
    // This is just a volume of air surrounding the detector so that
    // events (from calibration sources or cosmic rays) can be generated
    // on the outside.
    if (lab_walls_){
      // We want to simulate the walls (for muons in most cases).
      hallA_walls_->Construct();
      hallA_logic_ = hallA_walls_->GetLogicalVolume();
      G4double hallA_length = hallA_walls_->GetLSCHallALength();
      // Since the walls will be displaced need to make the
      // "lab" double sized to be sure.
      G4Box* lab_solid = new G4Box("LAB", hallA_length, hallA_length, hallA_length);
      G4Material *vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      lab_logic_ = new G4LogicalVolume(lab_solid, vacuum, "LAB");
      this->SetSpan(2 * hallA_length);
    }
    else {
      G4Box* lab_solid = new G4Box("LAB", lab_size_/2., lab_size_/2., lab_size_/2.);
      lab_logic_ = new G4LogicalVolume(lab_solid,
        G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "LAB");
    }
    lab_logic_->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Set this volume as the wrapper for the whole geometry
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(lab_logic_);

    // VESSEL (initialize first since it defines EL position)
    vessel_->SetELtoTPdistance(gate_tracking_plane_distance_);
    vessel_->Construct();
    G4LogicalVolume* vessel_logic = vessel_->GetLogicalVolume();
    G4LogicalVolume* vessel_internal_logic  = vessel_->GetInternalLogicalVolume();
    G4VPhysicalVolume* vessel_internal_phys = vessel_->GetInternalPhysicalVolume();
    G4ThreeVector vessel_displacement = shielding_->GetAirDisplacement(); // explained below
    gate_zpos_in_vessel_ = vessel_->GetELzCoord();

    // SHIELDING
    shielding_->Construct();
    shielding_->SetELzCoord(gate_zpos_in_vessel_);
    G4LogicalVolume* shielding_logic     = shielding_->GetLogicalVolume();
    G4LogicalVolume* shielding_air_logic = shielding_->GetAirLogicalVolume();

    // Recall that airbox is slighly displaced in Y dimension. In order to avoid
    // mistmatch with vertex generators, we place the vessel in the center of the world volume
    new G4PVPlacement(0, -vessel_displacement, vessel_logic,
                      "VESSEL", shielding_air_logic, false, 0);

    // INNER ELEMENTS
    inner_elements_->SetLogicalVolume(vessel_internal_logic);
    inner_elements_->SetPhysicalVolume(vessel_internal_phys);
    inner_elements_->SetELzCoord(gate_zpos_in_vessel_);
    inner_elements_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    inner_elements_->SetELtoTPdistance         (gate_tracking_plane_distance_);
    inner_elements_->Construct();

    // INNER COPPER SHIELDING
    ics_->SetLogicalVolume(vessel_internal_logic);
    ics_->SetELzCoord(gate_zpos_in_vessel_);
    ics_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    ics_->SetELtoTPdistance         (gate_tracking_plane_distance_);
    ics_->SetPortZpositions(vessel_->GetPortZpositions());
    ics_->Construct();

  //   // INNER TEFLON PANELS WITH FIBERS
  //   G4double vess_length = vessel_ -> GetLength();
  //   G4cout << "[Next100] Vessel length " << vess_length/1000 << " m " << G4endl;
   //
  //   // teflon panel /////////////////////////////////////////////////
  //   // materials
  //   G4Material* teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
  //   teflon->SetMaterialPropertiesTable(opticalprops::PTFE());
   //
  //   G4OpticalSurface* opsur_teflon =
  //     new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
  //   opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());
   //
  //   // dimensions
   //
  //   G4double panel_width_ = 170. * mm;
  //   G4double panel_thickness_ = 1. * mm;
   //
  //   G4Box* teflon_panel =
  //     new G4Box("TEFLON_PANEL", panel_width_/2., vess_length/2., panel_thickness_/2.);
  //   G4LogicalVolume* teflon_panel_logic =
  //     new G4LogicalVolume(teflon_panel, teflon, "TEFLON_PANEL");
   //
  //   new G4LogicalSkinSurface("TEFLON_PANEL_OPSURF", teflon_panel_logic, opsur_teflon);
   //
  //   teflon_panel_logic->SetVisAttributes(nexus::White());
  //   // teflon_panel_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
   //
  //   //
  //   // fiber ////////////////////////////////////////////////////
   //
  //   G4double fiber_diameter_ = 1. * mm;
  //   G4bool coated_ = true;
   //
  //   G4Material *this_fiber = nullptr;
  //   G4MaterialPropertiesTable *this_fiber_optical = nullptr;
  //   G4Material *this_coating = nullptr;
  //   G4MaterialPropertiesTable *this_coating_optical = nullptr;
   //
  //   if (fiber_type_ == "Y11") {
  //     this_fiber = materials::Y11();
  //     this_fiber_optical = opticalprops::Y11();
   //
  //     if (coated_) {
  //       this_coating = materials::TPB();
  //       this_coating_optical = opticalprops::TPB();
  //     }
   //
  //   } else if (fiber_type_ == "B2") {
   //
  //     this_fiber = materials::B2();
  //     this_fiber_optical = opticalprops::B2();
   //
  //     if (coated_) {
  //       this_coating = materials::TPH();
  //       this_coating_optical = opticalprops::TPH();
  //     }
   //
  //   } else {
  //     G4Exception("[FiberBarrel]", "Construct()",
  //                 FatalException, "Invalid fiber type, must be Y11 or B2");
  //   }
   //
   //
  //   fiber_ = new GenericWLSFiber(fiber_type_, true, fiber_diameter_, vess_length, true, coated_, this_coating, this_fiber, true);
   //
  //   fiber_->SetCoreOpticalProperties(this_fiber_optical);
  //   fiber_->SetCoatingOpticalProperties(this_coating_optical);
   //
  //   fiber_->Construct();
  //   G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
  //   if (fiber_type_ == "Y11")
  //     fiber_logic->SetVisAttributes(nexus::LightGreenAlpha());
  //   else if (fiber_type_ == "B2")
  //     fiber_logic->SetVisAttributes(nexus::LightBlueAlpha());
   //
  //   // GEOMETRY PARAMETERS /////////////////////////////////////////////
  //   G4double rot_angle;
   //
  //   G4double inner_rad = vessel_ -> GetInnerRadius();
  //   G4cout << "[Next100] Vessel inner radius " << inner_rad/10 << " cm radius" << G4endl;
   //
  //   // Teflon panels distance to the center
  //   G4double h = inner_rad - panel_thickness_/2.;
   //
  //   // Teflon panels angular separation
  //   G4double dif_theta = 2*std::atan(panel_width_/(2.*h));
   //
  //   G4int n_panels = floor(( 2 * M_PI) / dif_theta); // optimize the number of panels
   //
  //   G4int n_fibers = floor(panel_width_ / fiber_diameter_); // number of fibers per panel
  //   G4double dl_fib = panel_width_/n_fibers; // distance between fibers
   //
  //  G4cout << "[Next100] Using " << n_fibers << " fibers per panel"<< G4endl;
   //
  //  // Re-calculation of parameters for optimization
  //  dif_theta = ( 2 * M_PI) / n_panels; // re-calculate angular difference
  //  h = (panel_width_/2.)/(std::tan(dif_theta/2.)) + panel_thickness_/2. + fiber_diameter_; // re-calculate distance to the center
  //  G4cout << "[Next100] Using " << n_panels << " panels" << G4endl;
   //
  //   // Fibers distance to the center
  //  G4double hh = h - (fiber_diameter_/2. + panel_thickness_/2.);
   //
  //  // PLACEMENT /////////////////////////////////////////////
  //  // n_panels = 1;
   //
  //  for (G4int itheta=0; itheta < n_panels; itheta++) {
  //  // for (G4int itheta=0; itheta < 3; itheta++) {
   //
  //    // panels
  //    G4double theta = dif_theta * itheta;
  //    G4double x = h * std::cos(theta) * mm;
  //    G4double y = h * std::sin(theta) * mm;
  //    G4double phi = pi/2. + std::atan2(y, x);
  //    std::string label = std::to_string(itheta);
   //
  //    G4RotationMatrix* panel_rot = new G4RotationMatrix();
  //    rot_angle = pi/2.;
  //    panel_rot->rotateX(rot_angle);
  //    panel_rot->rotateY(phi);
  //    new G4PVPlacement(panel_rot, G4ThreeVector(x,y, 0.),
  //                      teflon_panel_logic, "PANEL-"+label, lab_logic_,
  //                      false, itheta, false);
   //
  //    // Relative positions of the fibers wrt the panel
  //    G4double x0_f = x*hh/h + (panel_width_/2. - fiber_diameter_/2.)*std::cos(phi);
  //    G4double y0_f = y*hh/h + (panel_width_/2. - fiber_diameter_/2.)*std::sin(phi);
   //
  //    for (G4int ii=0; ii < n_fibers; ii++) {
  //    // for (G4int ii=0; ii < 1; ii++) {
   //
  //        G4double xx_f = x0_f - dl_fib*ii*std::cos(phi);
  //        G4double yy_f = y0_f - dl_fib*ii*std::sin(phi);
   //
  //        std::string label2 = std::to_string(ii);
   //
  //        new G4PVPlacement(0, G4ThreeVector(xx_f, yy_f),
  //                          fiber_logic, "FIBER-" + label + label2, lab_logic_,
  //                          false, n_panels + ii, false);
  //      }
   //
  //  }



    //// PLACEMENT

    G4ThreeVector gate_pos(0., 0., -gate_zpos_in_vessel_);
    if (lab_walls_){
      G4ThreeVector castle_pos(0., hallA_walls_->GetLSCHallACastleY(),
                               hallA_walls_->GetLSCHallACastleZ());

      new G4PVPlacement(0, castle_pos, shielding_logic,
                        "LEAD_BOX", hallA_logic_, false, 0);
      new G4PVPlacement(0, gate_pos - castle_pos, hallA_logic_,
                        "Hall_A", lab_logic_, false, 0, false);
    }
    else {
      new G4PVPlacement(0, gate_pos, shielding_logic, "LEAD_BOX", lab_logic_, false, 0);
    }

    //// VERTEX GENERATORS
    lab_gen_ =
      new BoxPointSampler(lab_size_ - 1.*m, lab_size_ - 1.*m, lab_size_  - 1.*m, 1.*m,
                          G4ThreeVector(0., 0., 0.), 0);
  }


  G4ThreeVector Next100::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    // Air around shielding
    if (region == "LAB") {
      vertex = lab_gen_->GenerateVertex("INSIDE");
    }

    // Shielding regions
    else if ((region == "SHIELDING_LEAD")  ||
             (region == "SHIELDING_STEEL") ||
             (region == "INNER_AIR") ||
             (region == "EXTERNAL") ||
             (region == "SHIELDING_STRUCT") ||
             (region == "PEDESTAL") ||
             (region == "BUBBLE_SEAL") ||
             (region == "EDPM_SEAL")) {
      vertex = shielding_->GenerateVertex(region);
    }

    // Vessel regions
    else if ((region == "VESSEL")  ||
             (region == "PORT_1a") ||
             (region == "PORT_2a") ||
             (region == "PORT_1b") ||
             (region == "PORT_2b")) {
      vertex = vessel_->GenerateVertex(region);
    }

    // Inner copper shielding
    else if (region == "ICS"){
      vertex = ics_->GenerateVertex(region);
    }

    // Inner elements (photosensors' planes and field cage)
    else if ((region == "CENTER") ||
             (region == "ACTIVE") ||
             (region == "CATHODE_RING") ||
             (region == "BUFFER") ||
             (region == "XENON") ||
             (region == "LIGHT_TUBE") ||
             (region == "HDPE_TUBE") ||
             (region == "EL_GAP") ||
             (region == "EP_COPPER_PLATE") ||
             (region == "SAPPHIRE_WINDOW") ||
             (region == "OPTICAL_PAD") ||
             (region == "PMT_BODY") ||
             (region == "PMT") ||
             (region == "PMT_BASE") ||
             (region == "TP_COPPER_PLATE") ||
             (region == "SIPM_BOARD") ||
             (region == "DB_PLUG") ||
             (region == "EL_TABLE") ||
             (region == "FIELD_RING") ||
             (region == "GATE_RING") ||
             (region == "ANODE_RING") ||
             (region == "RING_HOLDER")) {
      vertex = inner_elements_->GenerateVertex(region);
    }

    else if (region == "AD_HOC") {
      // AD_HOC does not need to be shifted because it is passed by the user
      vertex = specific_vertex_;
      return vertex;
    }

    // Lab walls
    else if ((region == "HALLA_INNER") || (region == "HALLA_OUTER")){
      if (!lab_walls_)
        G4Exception("[Next100]", "GenerateVertex()", FatalException,
                    "This vertex generation region must be used with lab_walls == true!");
      vertex = hallA_walls_->GenerateVertex(region);
      while (vertex[1]<(-shielding_->GetHeight()/2.)){
        vertex = hallA_walls_->GenerateVertex(region);}
    }

    else {
      G4Exception("[Next100]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }

    G4ThreeVector displacement = G4ThreeVector(0., 0., -gate_zpos_in_vessel_);
    vertex = vertex + displacement;

    return vertex;
  }

} //end namespace nexus
