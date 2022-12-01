#ifndef FIBER_CYLINDER_H
#define FIBER_CYLINDER_H

#include "GeometryBase.h"
#include <G4MaterialPropertyVector.hh>
#include "GenericWLSFiber.h"

class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;
class CylinderPointSampler;

namespace nexus
{
    class CylinderPointSampler;
    class Fiber_cylinder: public GeometryBase
    {
        public:
        Fiber_cylinder();
        ~Fiber_cylinder();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        G4double radius_cyl_; //radius of the cylinder
        CylinderPointSampler* cyl_vertex_gen_; // this creates photons homogeneously in a cylinder
        G4GenericMessenger*   msg_;
    };
} // namespace nexus
#endif
