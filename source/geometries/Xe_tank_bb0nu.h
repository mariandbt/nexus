// ----------------------------------------------------------------------------
// nexus | Xe_tank_bb0nu.h
//
// Sphere filled with xenon.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef XE_SPHERE_H
#define XE_SPHERE_H

#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus {
class CylinderPointSampler2020;
}


namespace nexus {

  /// Cylindrical tank filled with xenon (liquid or gas)

  class Xe_tank_bb0nu: public GeometryBase
  {
  public:
    /// Constructor
    Xe_tank_bb0nu();
    /// Destructor
    ~Xe_tank_bb0nu();

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:
    G4bool liquid_;     ///< Whether xenon is liquid or not
    G4double pressure_; ///< Pressure (if gaseous state was selected)
    G4double radius_;   ///< Radius of the cylindric tank
    G4double length_;   ///< Length of the cylindric tank

    /// Vertexes random generator
    CylinderPointSampler2020* tank_vertex_gen_;

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
