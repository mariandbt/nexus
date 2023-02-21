#ifndef FIB_BOX_STRUCT_H
#define FIB_BOX_STRUCT_H

#include "GeometryBase.h"
#include <G4MaterialPropertyVector.hh>
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"

class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;

namespace nexus
{
    class Fib_box_struct: public GeometryBase
    {
        private:
        // fiber
        G4double diameter_;     // diameter of the cylindrical optical fibre
        G4double length_;     // length of the cylindrical optical fibre

        // box
        G4double box_xy_;     // outer side of the box structure
        G4double box_z_;     // outer length of the box structure
        G4double side_thickness_; // thickness of the box

        // sensor
        GenericPhotosensor* photo_sensor_;
        G4String sensor_type_; // SiPM, PMT, PERFECT, ...
        G4double sensor_width_;   // Width of rectangular fiber sensors
        G4double sensor_height_;   // Height of rectangular fiber sensors
        G4double sensor_thickness_;   // (Thickness set to a fix value of 1 mm)
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
