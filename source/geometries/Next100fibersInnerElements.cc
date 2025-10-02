// ----------------------------------------------------------------------------
// nexus | Next100fibersInnerElements.cc
//
// Inner elements of the NEXT-100 detector. They include the field cage,
// the tracking plane.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100fibersInnerElements.h"
#include "Next100fibersFieldCage.h"
#include "Next100fibersTrackingPlane.h"

#include <G4GenericMessenger.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Material.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace CLHEP;


namespace nexus {


  Next100fibersInnerElements::Next100fibersInnerElements(G4double grid_thickn):
    GeometryBase(),
    mother_logic_(nullptr),
    mother_phys_ (nullptr),
    gas_(nullptr),
    field_cage_    (new Next100fibersFieldCage(grid_thickn)),
    tracking_plane_(new Next100fibersTrackingPlane()),
    msg_(nullptr)
  {
    // Messenger
    msg_ = new G4GenericMessenger(this, "/Geometry/Next100fibers/",
                                  "Control commands of geometry Next100fibers.");
  }


  void Next100fibersInnerElements::SetLogicalVolume(G4LogicalVolume* mother_logic)
  {
    mother_logic_ = mother_logic;
  }


  void Next100fibersInnerElements::SetPhysicalVolume(G4VPhysicalVolume* mother_phys)
  {
    mother_phys_ = mother_phys;
  }


  void Next100fibersInnerElements::Construct()
  {
    G4ThreeVector coord_origin = GetCoordOrigin();
    // Reading mother material
    gas_ = mother_logic_->GetMaterial();
    pressure_ =    gas_->GetPressure();
    temperature_ = gas_->GetTemperature();

    // Field Cage
    field_cage_->SetMotherLogicalVolume(mother_logic_);
    field_cage_->SetMotherPhysicalVolume(mother_phys_);
    field_cage_->SetCoordOrigin(coord_origin);
    field_cage_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    field_cage_->SetSiPMPitch(tracking_plane_->GetSiPMPitch());
    field_cage_->Construct();

    // Tracking plane
    tracking_plane_->SetMotherPhysicalVolume(mother_phys_);
    tracking_plane_->SetCoordOrigin(coord_origin);
    tracking_plane_->SetELtoTPdistance(gate_tracking_plane_distance_);
    tracking_plane_->Construct();

    tracking_plane_->GetSiPMPosInGas(sipm_pos_);
  }


  Next100fibersInnerElements::~Next100fibersInnerElements()
  {
    delete field_cage_;
    delete tracking_plane_;
  }


  G4ThreeVector Next100fibersInnerElements::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    // Field Cage regions
    if ((region == "CENTER") ||
        (region == "ACTIVE") ||
        (region == "CATHODE_RING") ||
        (region == "CATHODE_SURF") ||
        (region == "BUFFER") ||
        (region == "XENON")  ||
        (region == "S2_PMT_LT") ||
        (region == "S2_SIPM_PSF") ||
        (region == "LIGHT_TUBE") ||
        (region == "HDPE_TUBE") ||
        (region == "FIELD_RING") ||
        (region == "GATE_RING") ||
        (region == "ANODE_RING") ||
        (region == "RING_HOLDER")) {
      vertex = field_cage_->GenerateVertex(region);
    }
    // Teflon Plane regions
    else if ((region == "EP_COPPER_PLATE") ||
             (region == "SAPPHIRE_WINDOW") ||
             (region == "OPTICAL_PAD") ||
             (region == "PMT") ||
             (region == "PMT_BODY") ||
             (region == "PMT_BASE")) {
      // vertex = energy_plane_->GenerateVertex(region);
    }
    // Tracking Plane regions
    else if ((region == "TP_COPPER_PLATE") ||
             (region == "SIPM_BOARD") ||
             (region == "DB_PLUG")) {
      vertex = tracking_plane_->GenerateVertex(region);
    }
    else {
      G4Exception("[Next100fibersInnerElements]", "GenerateVertex()", FatalException,
        "Unknown vertex generation region!");
    }

    return vertex;
  }

} // end namespace nexus
