// ----------------------------------------------------------------------------
// nexus | NextHDDEMOOpticalGeometry.h
//
// This class builds a simplified version of the NEXT-100 geometry, where
// only the inner elements are instantiated.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_OPT_GEO_H
#define NEXTHDDEMO_OPT_GEO_H

#include <G4ThreeVector.hh>
#include "GeometryBase.h"

class G4GenericMessenger;


namespace nexus {

  class NextHDDEMOInnerElements;

  class NextHDDEMOOpticalGeometry : public GeometryBase
  {

  public:
    ///Constructor
    NextHDDEMOOpticalGeometry();

    /// Destructor
    ~NextHDDEMOOpticalGeometry();

    /// Returns a vertex in a region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    /// Builder
    void Construct();


  private:

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    G4double gate_tracking_plane_distance_, gate_sapphire_wdw_distance_;
    G4double pressure_;
    G4double temperature_;
    G4double sc_yield_;
    G4double e_lifetime_;

    // Vertex decided by user
    G4ThreeVector specific_vertex_;

    G4String gas_;

    NextHDDEMOInnerElements* inner_elements_;

    // Relative position of the gate in its mother volume
    G4double gate_zpos_in_gas_;

  };

} // end namespace nexus

#endif
