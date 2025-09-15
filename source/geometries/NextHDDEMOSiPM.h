// -----------------------------------------------------------------------------
//  nexus | NextHDDEMOSiPM.h
//
//  Geometry of the Hamamatsu MPPC S13372-1350TE, the model of
//  silicon photomultiplier (SiPM) used in the NEXT-100 detector.
//
//  The NEXT Collaboration
// -----------------------------------------------------------------------------

#ifndef NEXTHDDEMO_SIPM_H
#define NEXTHDDEMO_SIPM_H

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

namespace nexus {

  class NextHDDEMOSiPM: public GeometryBase
  {
  public:
    // Constructor
    NextHDDEMOSiPM();
    /// Destructor
    ~NextHDDEMOSiPM();

    // Return dimensions of the SiPM
    G4ThreeVector GetDimensions() const;

    // Invoke this method to build the volumes of the geometry
    void Construct() override;

    // Needed settings for correct numbering
    void SetSensorDepth         (G4int sensor_depth);
    void SetMotherDepth         (G4int mother_depth);
    void SetNamingOrder         (G4int naming_order);
    void SetTimeBinning         (G4double time_binning);
    void SetSiPMCoatingThickness(G4double coating_thickn);
    void SetVisibility          (G4bool visibility);

  private:
    G4ThreeVector dimensions_;

    G4int    sensor_depth_;
    G4int    mother_depth_;
    G4int    naming_order_;
    G4double time_binning_;

    G4double coating_thickn_;
    G4bool visibility_;

  };

  inline void NextHDDEMOSiPM::SetTimeBinning(G4double time_binning)
  { time_binning_ = time_binning; }

  inline void NextHDDEMOSiPM::SetSensorDepth(G4int sensor_depth)
  { sensor_depth_ = sensor_depth; }

  inline void NextHDDEMOSiPM::SetMotherDepth(G4int mother_depth)
  { mother_depth_ = mother_depth; }

  inline void NextHDDEMOSiPM::SetNamingOrder(G4int naming_order)
  { naming_order_ = naming_order; }

  inline void NextHDDEMOSiPM::SetSiPMCoatingThickness(G4double coating_thickn)
  { coating_thickn_ = coating_thickn; }

  inline void NextHDDEMOSiPM::SetVisibility(G4bool visibility)
  { visibility_ = visibility; }

}

#endif
