#ifndef FIBER_MORE_H
#define FIBER_MORE_H

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
    class Fiber_more: public GeometryBase
    {
        public:
        Fiber_more();
        ~Fiber_more();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        // CylinderPointSampler* cyl_vertex_gen_; // this creates photons homogeneously in a cylinder, so it is used for the tests but for the simulation we want to create the phton outside the fiber
        G4GenericMessenger*   msg_;
    };
} // namespace nexus
#endif
