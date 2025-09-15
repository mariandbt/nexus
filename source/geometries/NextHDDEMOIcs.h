// ----------------------------------------------------------------------------
// nexus | NextHDDEMOIcs.h
//
// Inner copper shielding of the NEXT-100 detector.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_ICS_H
#define NEXTHDDEMO_ICS_H

#include "GeometryBase.h"

#include <G4Navigator.hh>

class G4GenericMessenger;


namespace nexus {

  class CylinderPointSampler2020;

  class NextHDDEMOIcs: public GeometryBase
  {
  public:
    /// Constructor
    NextHDDEMOIcs();

    /// Destructor
    ~NextHDDEMOIcs();

    /// Sets the Logical Volume where ICS will be placed
    void SetLogicalVolume(G4LogicalVolume* mother_logic);

    void SetELtoTPdistance(G4double);
    void SetELtoSapphireWDWdistance(G4double);
    void SetPortZpositions(G4double port_positions[]);

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    /// Builder
    void Construct();

  private:
    // Mother Logical Volume of the ICS
    G4LogicalVolume* mother_logic_;

    // Dimensions
    G4double gate_tp_distance_, gate_sapphire_wdw_dist_;
    G4double in_rad_, thickness_;
    G4double port_z_1a_, port_z_2a_, port_z_1b_, port_z_2b_;

    // Visibility of the shielding
    G4bool visibility_;

    // Vertex generator
    CylinderPointSampler2020* ics_gen_;

    // Geometry Navigator
    G4Navigator* geom_navigator_;

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

  };

  inline void NextHDDEMOIcs::SetELtoTPdistance(G4double distance){
    gate_tp_distance_ = distance;
  }

  inline void NextHDDEMOIcs::SetELtoSapphireWDWdistance(G4double distance){
    gate_sapphire_wdw_dist_ = distance;
  }

} // end namespace nexus

#endif
