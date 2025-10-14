// ----------------------------------------------------------------------------
// nexus | NextHDDEMOInnerElements.cc
//
// Inner elements of the NEXT-HDDemo detector. They include the field cage,
// the energy and the tracking plane.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "NextHDDEMOInnerElements.h"
#include "NextHDDEMOFieldCage.h"
#include "NextHDDEMOTrackingPlane.h"

#include <G4GenericMessenger.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Material.hh>

#include <CLHEP/Units/SystemOfUnits.h>

using namespace CLHEP;


namespace nexus {


  NextHDDEMOInnerElements::NextHDDEMOInnerElements():
    GeometryBase(),
    mother_logic_(nullptr),
    mother_phys_ (nullptr),
    gas_(nullptr),
    field_cage_    (new NextHDDEMOFieldCage()),
    tracking_plane_(new NextHDDEMOTrackingPlane()),
    msg_(nullptr)
  {
    // Messenger
    msg_ = new G4GenericMessenger(this, "/Geometry/NextHDDEMO/",
                                  "Control commands of geometry NextHDDEMO.");
  }


  void NextHDDEMOInnerElements::SetLogicalVolume(G4LogicalVolume* mother_logic)
  {
    mother_logic_ = mother_logic;
  }


  void NextHDDEMOInnerElements::SetPhysicalVolume(G4VPhysicalVolume* mother_phys)
  {
    mother_phys_ = mother_phys;
  }


  void NextHDDEMOInnerElements::Construct()
  {
    G4cout << "[NextHDDEMOInnerElements] Constructor called" << G4endl;
    // Position in Z of the beginning of the drift region
    G4double gate_zpos = GetELzCoord();
    // Reading mother material
    gas_ = mother_logic_->GetMaterial();
    pressure_ =    gas_->GetPressure();
    temperature_ = gas_->GetTemperature();

    // Field Cage
    field_cage_->SetMotherLogicalVolume(mother_logic_);
    field_cage_->SetMotherPhysicalVolume(mother_phys_);
    field_cage_->SetELzCoord(gate_zpos);
    field_cage_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    field_cage_->Construct();

    // Tracking plane
    tracking_plane_->SetMotherPhysicalVolume(mother_phys_);
    tracking_plane_->SetELzCoord(gate_zpos);
    tracking_plane_->SetELtoTPdistance(gate_tracking_plane_distance_);
    tracking_plane_->Construct();
  }


  NextHDDEMOInnerElements::~NextHDDEMOInnerElements()
  {
    delete field_cage_;
    // delete energy_plane_;
    delete tracking_plane_;
  }


  G4ThreeVector NextHDDEMOInnerElements::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    // Field Cage regions
    if ((region == "CENTER") ||
    (region == "SECTION_AREA") ||
    (region == "SEGMENT") ||
    (region == "SECTOR_AREA") ||
    (region == "SECTOR_VOL") ||
        (region == "ACTIVE") ||
        (region == "CATHODE_RING") ||
        (region == "BUFFER") ||
        (region == "XENON")  ||
        (region == "EL_GAP") ||
        (region == "LIGHT_TUBE") ||
        (region == "HDPE_TUBE") ||
        (region == "FIELD_RING") ||
        (region == "GATE_RING") ||
        (region == "ANODE_RING") ||
        (region == "RING_HOLDER")) {
      vertex = field_cage_->GenerateVertex(region);
    }
    // Tracking Plane regions
    else if ((region == "TP_COPPER_PLATE") ||
             (region == "SIPM_BOARD") ||
             (region == "DB_PLUG")) {
      vertex = tracking_plane_->GenerateVertex(region);
    }
    else {
      G4Exception("[NextHDDEMOInnerElements]", "GenerateVertex()", FatalException,
        "Unknown vertex generation region!");
    }

    return vertex;
  }

} // end namespace nexus
