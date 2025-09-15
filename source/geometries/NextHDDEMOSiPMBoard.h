// -----------------------------------------------------------------------------
// nexus | NextHDDEMOSiPMBoard.h
//
// Geometry of the NEXT-100 SiPM board, consisting of an 8x8 array of
// silicon photomultipliers (1.3x1.3 mm2 of active area) mounted on a Kapton
// board covered with a TPB-coated teflon mask.
//
// The NEXT Collaboration
// -----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_SIPM_BOARD_H
#define NEXTHDDEMO_SIPM_BOARD_H

#include "GeometryBase.h"
#include <G4ThreeVector.hh>
#include <vector>

class G4VPhysicalVolume;
class G4GenericMessenger;

namespace nexus {

  class BoxPointSampler;
  class NextHDDEMOSiPM;

  // Geometry of the 8x8 SiPM boards used in the tracking plane of NEXT-100

  class NextHDDEMOSiPMBoard: public GeometryBase
  {
  public:
    // Default constructor
    NextHDDEMOSiPMBoard();
    // Destructor
    ~NextHDDEMOSiPMBoard();
    //
    void SetMotherPhysicalVolume(G4VPhysicalVolume*);
    //
    void Construct() override;
    //
    G4ThreeVector GenerateVertex(const G4String&) const override;

    G4double GetSize() const;
    G4double GetThickness() const;

    const std::vector<G4ThreeVector>& GetSiPMPositions() const;

  private:
    G4GenericMessenger* msg_;
    G4double size_, pitch_, margin_;
    G4double board_thickness_, mask_thickness_;
    G4double time_binning_;
    std::vector<G4ThreeVector> sipm_positions_;
    G4bool   visibility_, sipm_visibility_;
    G4VPhysicalVolume*  mpv_;
    BoxPointSampler*    vtxgen_;
    NextHDDEMOSiPM* sipm_;
  };

  inline void NextHDDEMOSiPMBoard::SetMotherPhysicalVolume(G4VPhysicalVolume* p)
  { mpv_ = p;}

  inline G4double NextHDDEMOSiPMBoard::GetSize() const
  { return size_; }

  inline G4double NextHDDEMOSiPMBoard::GetThickness() const
  { return (board_thickness_ + mask_thickness_); }

  inline const std::vector<G4ThreeVector>& NextHDDEMOSiPMBoard::GetSiPMPositions() const
  { return sipm_positions_; }

} // namespace nexus

#endif
