#ifndef FIB_BOX_STRUCT_H
#define FIB_BOX_STRUCT_H

#include "GeometryBase.h"
#include <G4MaterialPropertyVector.hh>
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"
// #include "SiPM11_eff.h"

class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;

namespace nexus
{
    class Fib_box_struct: public GeometryBase
    {
        private:
        // fiber
        G4double radius_;     // radius of the cylindrical optical fibre
        G4double length_;     // length of the cylindrical optical fibre

        // box
        G4double box_xy_;     // outer side of the box structure
        G4double box_z_;     // outer length of the box structure
        G4double side_thickness; // thickness of the box

        // sensor
        // SiPM11_eff* sipm_;   // SiPM
        GenericPhotosensor* photo_sensor_;
        G4bool sensor_visibility_;

        G4GenericMessenger*   msg_;

        public:
        Fib_box_struct();
        ~Fib_box_struct();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;
    };
} // namespace nexus
#endif
