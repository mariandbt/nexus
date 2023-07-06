// ----------------------------------------------------------------------------
// nexus | FibBarrMeth.cc
//
// Cylinder containing optical fibers with a methacrylate window
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------


#ifndef FibBarrMeth_H
#define FibBarrMeth_H

#include "GeometryBase.h"
#include "GenericWLSFiber.h"
#include "GenericCircularPhotosensor.h"
#include "CylinderPointSampler2020.h"
#include "MaterialsList.h"

#include <G4MaterialPropertyVector.hh>


class G4Material;
class G4GenericMessenger;
class G4MaterialPropertiesTable;
class CylinderPointSampler2020;

namespace nexus
{
    class CylinderPointSampler2020;
    class FibBarrMeth: public GeometryBase
    {
        public:
        FibBarrMeth();
        ~FibBarrMeth();
        void Construct();
        G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
        // fiber
        G4String fiber_type_;   // Y11 or B2
        G4bool coated_;         // true or false for the fibers to be coated
        G4double diameter_;     // diameter of the cylindrical optical fibre
        G4double length_;       // length of the cylindrical optical fibre

        // cylinder
        G4double radius_cyl_;         //radius of the cylinder

        // // sensor
        // GenericCircularPhotosensor* photo_sensor_;
        // G4bool sensor_visibility_;

        CylinderPointSampler2020* cyl_vertex_gen_; // this creates photons homogeneously in a cylinder
        G4GenericMessenger*   msg_;


    };
} // namespace nexus
#endif
