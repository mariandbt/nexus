// ----------------------------------------------------------------------------
// nexus | Fib_pan_meth_MEDIUM.h
//
// Box-shaped box of material with a coating.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef Fib_pan_meth_MEDIUM_H
#define Fib_pan_meth_MEDIUM_H

#include "GeometryBase.h"
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"
#include "CylinderPointSampler2020.h"
#include "MaterialsList.h"

class G4Material;
class G4GenericMessenger;
namespace nexus { class SpherePointSampler; }


namespace nexus {

  /// Fiber Barrel

  class Fib_pan_meth_MEDIUM: public GeometryBase
  {
  public:
    /// Constructor
    Fib_pan_meth_MEDIUM();
    /// Destructor
    ~Fib_pan_meth_MEDIUM();

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:
    G4double world_z_;             // World dimensions
    G4double world_xy_;

    G4bool liquid_;               // Whether xenon is liquid or not
    G4double pressure_;           // Pressure (if gaseous state was selected)

    G4double radius_;             // Cylinder radius
    G4double fiber_diameter_;
    G4double panel_length_;             // Cylinder length
    G4String fiber_type_;
    G4bool coated_;

    // cylinder
    G4double teflon_thickness_;   // thickness of the outer teflon cover
    G4bool caps_visibility_;
    G4bool teflon_visibility_;

    G4double panel_width_;

    // methacrylate
    G4bool methacrylate_;
    G4double window_thickness_;

    // sensor
    GenericPhotosensor* photo_sensor_;
    G4String sensor_type_;        // SiPM, PMT, PERFECT, ...
    G4bool sensor_visibility_;

    GenericWLSFiber* fiber_;
    CylinderPointSampler2020* cyl_vertex_gen_;
    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
