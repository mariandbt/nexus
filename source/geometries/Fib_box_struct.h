#ifndef FIB_BOX_STRUCT_H
#define FIB_BOX_STRUCT_H

#include "GeometryBase.h"
#include <G4MaterialPropertyVector.hh>
#include "GenericWLSFiber.h"
#include "SiPM11_eff.h"

class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;

namespace nexus
{
    class Fib_box_struct: public GeometryBase
    {
        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        G4double box_xy_;     //side of the box structure
        G4double box_z_;     //length of the box structure
        SiPM11_eff* sipm_;   // SiPM
        G4GenericMessenger*   msg_;

        public:
        Fib_box_struct();
        ~Fib_box_struct();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;
    };
} // namespace nexus
#endif
