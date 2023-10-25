// ----------------------------------------------------------------------------
// nexus | Fiber_panel.h
//
// Box-shaped box of material with a coating.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef Fiber_panel_H
#define Fiber_panel_H

#include "GeometryBase.h"
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"
#include "MaterialsList.h"

class G4Material;
class G4GenericMessenger;


namespace nexus {

  /// Fiber Panel

  class Fiber_panel: public GeometryBase
  {
  public:
    /// Constructor
    Fiber_panel(
      G4String name,
      G4double panel_x_, G4double panel_y_, G4double panel_z_,
      G4String fiber_type_, G4double fiber_diameter_, G4bool coated_,
      G4String sens_type_, G4int n_sens_
    );
    /// Destructor
    ~Fiber_panel();

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

    //
    G4double GetWidth()       const;
    G4double GetHeight()      const;
    G4double GetThickness()   const;
    const G4String& GetName() const;

    void SetVisibility           (G4bool visibility);
    void SetOpticalProperties    (G4MaterialPropertiesTable* mpt);
    void SetTimeBinning          (G4double time_binning);
    void SetSensorDepth          (G4int sensor_depth);
    void SetMotherDepth          (G4int mother_depth);
    void SetNamingOrder          (G4int naming_order);


  private:

    G4String name_;

    /// panel
    G4double panel_x_;             // panel width
    G4double panel_y_;             // panel length
    G4double panel_z_;             // panel thickness

    /// aluminized fiber
    GenericWLSFiber* fiber_;
    G4String fiber_type_;         // Y11 or B2
    G4double fiber_diameter_;     // diameter of the round fiber
    G4bool coated_;               // whether or not to coat the fibers
    G4double fiber_end_z_;        // thickness of the aluminium end of the fiber

    // sensor
    GenericPhotosensor* photo_sensor_;
    G4String sens_type_;          // SiPM, PMT, PERFECT, ...
    G4String sens_z_;             // sensor thinckness
    G4int n_sens_;                // number of sensor on the panel

    // visibility
    G4bool visibility_;           // make geometry (not) visible

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
