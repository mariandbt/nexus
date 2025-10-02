// ----------------------------------------------------------------------------
// nexus | Next100fibersInnerElements.h
//
// Inner elements of the NEXT-100 detector. They include the field cage,
// and the tracking plane.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXT100_INNER_ELEMENTS_H
#define NEXT100_INNER_ELEMENTS_H

#include "GeometryBase.h"

#include <G4ThreeVector.hh>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4GenericMessenger;


namespace nexus {

  class Next100fibersFieldCage;
  class Next100fibersTrackingPlane;

  class Next100fibersInnerElements : public GeometryBase
  {

  public:
    ///Constructor
    Next100fibersInnerElements(G4double grid_thickn);

    /// Destructor
    ~Next100fibersInnerElements();

    /// Set the logical and physical volume that encloses the entire geometry
    void SetLogicalVolume(G4LogicalVolume*);
    void SetPhysicalVolume(G4VPhysicalVolume*);
    void SetELtoTPdistance(G4double);
    void SetELtoSapphireWDWdistance(G4double);

    /// Return the relative position respect to the rest of NEXT100 geometry
    G4ThreeVector GetPosition() const;

    /// Return the positions of the SiPMs in their mother volume (gas)
    std::vector<G4ThreeVector> GetSiPMPosInGas() const;

    /// Return the positions of the PMTs in their mother volume (gas)
    std::vector<G4ThreeVector> GetPMTPosInGas() const;

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
    Next100fibersFieldCage*     field_cage_;
    Next100fibersTrackingPlane* tracking_plane_;

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    // Positions of the SiPMs in their mother volume (gas)
    std::vector<G4ThreeVector> sipm_pos_;

    // Positions of the PMTs in their mother volume (gas)
    std::vector<G4ThreeVector> pmt_pos_;

  };

  inline void Next100fibersInnerElements::SetELtoTPdistance(G4double distance){
    gate_tracking_plane_distance_ = distance;
  }

  inline void Next100fibersInnerElements::SetELtoSapphireWDWdistance(G4double distance){
    gate_sapphire_wdw_distance_ = distance;
  }

  inline std::vector<G4ThreeVector> Next100fibersInnerElements::GetSiPMPosInGas() const
  {
    return sipm_pos_;
  }

 inline std::vector<G4ThreeVector> Next100fibersInnerElements::GetPMTPosInGas() const
  {
    return pmt_pos_;
  }

} // end namespace nexus

#endif
