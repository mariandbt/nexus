// ----------------------------------------------------------------------------
// nexus | noOpticalTrackingAction.h
//
// This class stores in memory the trajectories of all particles, except optical photons,
// including ionization electrons, with the relevant tracking information that will be
// saved to the output file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NO_OPTICAL_TRACKING_ACTION_H
#define NO_OPTICAL_TRACKING_ACTION_H

#include <G4UserTrackingAction.hh>

class G4Track;


namespace nexus {

  // General-purpose user tracking action

  class noOpticalTrackingAction: public G4UserTrackingAction
  {
  public:
    /// Constructor
    noOpticalTrackingAction();
    /// Destructor
    virtual ~noOpticalTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
  };

}

#endif
