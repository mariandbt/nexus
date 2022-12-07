// ----------------------------------------------------------------------------
// nexus | SiPM11_eff.h
//
// Geometry of a Hamamatsu 1x1 mm2 SiPM.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef SILICON_PM_11_EFF_H
#define SILICON_PM_11_EFF_H

#include "GeometryBase.h"
#include <G4ThreeVector.hh>

class G4GenericMessenger;

namespace nexus {


  /// Geometry of the Hamamatsu surface-mounted 1x1 mm2 MPPC (SiPM)

  class SiPM11_eff: public GeometryBase
  {
  public:
    /// Constructor
    SiPM11_eff();
    /// Destructor
    ~SiPM11_eff();

    /// Return dimensions of the SiPM
    G4ThreeVector GetDimensions() const;

    /// Invoke this method to build the volumes of the geometry
    void Construct();

  private:
    G4ThreeVector dimensions_; ///< external dimensions of the SiPM11_eff

    // Visibility of the tracking plane
    G4bool sipm_visibility_;

     // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

  };


} // end namespace nexus

#endif
