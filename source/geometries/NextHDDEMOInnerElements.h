// ----------------------------------------------------------------------------
// nexus | NextHDDEMOInnerElements.h
//
// Inner elements of the NEXT-HDDemo detector. They include the field cage,
// the energy and the tracking plane.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_INNER_ELEMENTS_H
#define NEXTHDDEMO_INNER_ELEMENTS_H

#include <G4ThreeVector.hh>
#include <vector>
#include "GeometryBase.h"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4GenericMessenger;


namespace nexus {

  class NextHDDEMOFieldCage;
  class NextHDDEMOTrackingPlane;

  class NextHDDEMOInnerElements : public GeometryBase
  {

  public:
    ///Constructor
    NextHDDEMOInnerElements();

    /// Destructor
    ~NextHDDEMOInnerElements();

    /// Set the logical and physical volume that encloses the entire geometry
    void SetLogicalVolume(G4LogicalVolume*);
    void SetPhysicalVolume(G4VPhysicalVolume*);
    void SetELtoTPdistance(G4double);
    void SetELtoSapphireWDWdistance(G4double);

    /// Return the relative position respect to the rest of NEXTHDDEMO geometry
    G4ThreeVector GetPosition() const;

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;

    /// Builder
    void Construct();


  private:

    G4double gate_sapphire_wdw_distance_;
    G4double gate_tracking_plane_distance_;


    G4LogicalVolume* mother_logic_;
    G4VPhysicalVolume* mother_phys_;
    G4Material* gas_;

    G4double pressure_;
    G4double temperature_;

    // Detector parts
    NextHDDEMOFieldCage*     field_cage_;
    NextHDDEMOTrackingPlane* tracking_plane_;

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

  };

  inline void NextHDDEMOInnerElements::SetELtoTPdistance(G4double distance){
    gate_tracking_plane_distance_ = distance;
  }

  inline void NextHDDEMOInnerElements::SetELtoSapphireWDWdistance(G4double distance){
    gate_sapphire_wdw_distance_ = distance;
  }

} // end namespace nexus

#endif
