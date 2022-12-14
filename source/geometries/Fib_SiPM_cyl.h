#ifndef FIB_SIPM_CYL_H
#define FIB_SIPM_CYL_H

#include "GeometryBase.h"
#include <G4MaterialPropertyVector.hh>
#include "GenericWLSFiber.h"
#include "SiPM11_eff.h"

class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;
class CylinderPointSampler;

namespace nexus
{
    class CylinderPointSampler;
    class Fib_SiPM_cyl: public GeometryBase
    {
        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        G4double radius_cyl_; //radius of the cylinder
        SiPM11_eff* sipm_;   // SiPM
        CylinderPointSampler* cyl_vertex_gen_; // this creates photons homogeneously in a cylinder
        G4GenericMessenger*   msg_;

        public:
        Fib_SiPM_cyl();
        ~Fib_SiPM_cyl();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;
    };
} // namespace nexus
#endif
