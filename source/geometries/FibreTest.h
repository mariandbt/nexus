#ifndef FIBRE_TEST_H
#define FIBRE_TEST_H

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
    class FibreTest: public GeometryBase
    {
        public:
        FibreTest();
        ~FibreTest();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        CylinderPointSampler* cyl_vertex_gen_;
        G4GenericMessenger*   msg_;
    };
} // namespace nexus
#endif
