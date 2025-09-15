// ----------------------------------------------------------------------------
// nexus | NextHDDEMO.h
//
// Main class that constructs the geometry of the NEXT-100 detector.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_H
#define NEXTHDDEMO_H

#include "GeometryBase.h"

// // Marian's adenda
// #include "GenericWLSFiber.h"
// #include "MaterialsList.h"
// //

class G4LogicalVolume;
class G4GenericMessenger;

namespace nexus {class BoxPointSampler;}


namespace nexus {

  class BoxPointSampler;
  class NextHDDEMOShielding;
  class NextHDDEMOVessel;
  class NextHDDEMOIcs;
  class NextHDDEMOInnerElements;
  class LSCHallA;

  class NextHDDEMO: public GeometryBase
  {
  public:
    /// Constructor
    NextHDDEMO();

    /// Destructor
    ~NextHDDEMO();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;


  private:
    void BuildLab();
    void Construct();


  private:
    // Detector dimensions
    const G4double lab_size_;          /// Size of the air box containing the detector
    const G4double gate_tracking_plane_distance_, gate_sapphire_wdw_distance_;

    // Pointers to logical volumes
    G4LogicalVolume* lab_logic_;
    G4LogicalVolume* buffer_gas_logic_;
    G4LogicalVolume* hallA_logic_;

    // Detector parts
    LSCHallA* hallA_walls_;
    NextHDDEMOShielding* shielding_;
    NextHDDEMOVessel*    vessel_;
    NextHDDEMOIcs*       ics_;
    NextHDDEMOInnerElements* inner_elements_;

    BoxPointSampler* lab_gen_; ///< Vertex generator

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    /// Specific vertex for AD_HOC region
    G4ThreeVector specific_vertex_;

    /// Position of gate in its mother volume
    G4double gate_zpos_in_vessel_;

    /// Whether or not to build LSC HallA.
    G4bool lab_walls_;

    // // WSL fibers (Y11 or B2)
    // GenericWLSFiber* fiber_;
    // G4String fiber_type_;
  };

} // end namespace nexus

#endif
