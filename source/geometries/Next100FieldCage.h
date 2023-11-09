// ----------------------------------------------------------------------------
// nexus | Next100FieldCage.h
//
// Geometry of the NEXT-100 field cage.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXT100_FIELDCAGE_H
#define NEXT100_FIELDCAGE_H

#include "GeometryBase.h"
#include <vector>

// Marian's adenda
#include "GenericWLSFiber.h"
#include "GenericPhotosensor.h"
#include "MaterialsList.h"
//

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GenericMessenger;
class G4Navigator;

namespace nexus {

  class CylinderPointSampler2020;
  class SegmentPointSampler;


  class Next100FieldCage: public GeometryBase
  {
  public:
    Next100FieldCage();
    ~Next100FieldCage();
    void Construct() override;
    G4ThreeVector GenerateVertex(const G4String& region) const override;

    G4ThreeVector GetActivePosition() const;

    void SetMotherLogicalVolume(G4LogicalVolume* mother_logic);
    void SetMotherPhysicalVolume(G4VPhysicalVolume* mother_phys);
    void SetELtoSapphireWDWdistance(G4double);

  private:
    void DefineMaterials();
    void BuildActive();
    void BuildCathode();
    void BuildBuffer();
    void BuildELRegion();
    void BuildFiberBarrel();
    void BuildLightTube();
    void BuildFieldCage();

    // Dimensions
    G4double gate_sapphire_wdw_dist_;
    const G4double active_diam_;
    const G4double cathode_int_diam_, cathode_ext_diam_, cathode_thickn_, grid_thickn_;
    const G4double teflon_drift_length_, teflon_total_length_, teflon_thickn_;
    // const G4int n_panels_;
    const G4double el_gap_length_;
    const G4double gate_teflon_dist_, gate_ext_diam_, gate_int_diam_, gate_ring_thickn_;
    const G4double hdpe_tube_int_diam_, hdpe_tube_ext_diam_, hdpe_length_;
    const G4double ring_ext_diam_, ring_int_diam_, ring_thickn_, drift_ring_dist_, buffer_ring_dist_;
    const G4double holder_x_, holder_long_y_, holder_short_y_;
    const G4double tpb_thickn_;
    const G4double overlap_;

    // Fiber Barrel
    //// fibers
    GenericWLSFiber* fiber_; // WSL fibers (Y11 or B2)
    G4String fiber_type_; // WSL fibers (Y11 or B2)
    G4double fiber_diameter_;
    G4double fiber_end_z;
    G4double fiber_length;
    G4bool coated_; // wheter the fibers have coating or not
    //// sensor
    GenericPhotosensor* photo_sensor_;
    G4String sensor_type_;        // SiPM, PMT, PERFECT, ...
    G4bool sensor_visibility_;
    G4double sens_z;
    //// barrel
    G4bool cap_visibility_;
    G4bool panels_visibility_;
    G4double panel_width_;
    G4double panel_thickness_;
    G4double panel_length_;
    //// relative z-positions to the panels (reference z position)
    G4double z, z_f, z_fend, z_s;
    //// relative positions to the center
    G4double h, hh;
    //// space between elements
    G4double dif_theta, dl_fib, dl_sens;
    //// number of elements
    G4int n_panels, n_fibers, n_sensors;


    // Diffusion constants
    G4double drift_transv_diff_, drift_long_diff_;
    G4double ELtransv_diff_; ///< transversal diffusion in the EL gap
    G4double ELlong_diff_; ///< longitudinal diffusion in the EL gap
    // Electric field
    G4bool elfield_;
    G4double ELelectric_field_; ///< electric field in the EL region
    // Transparencies of grids
    const G4double cath_grid_transparency_, el_grid_transparency_;

    //Step size
    G4double max_step_size_;

    // Visibility of the geometry
    G4bool visibility_;
    // Verbosity of the geometry
    G4bool verbosity_;

    G4double active_length_, buffer_length_;
    G4double teflon_buffer_length_;
    G4double teflon_drift_zpos_,teflon_buffer_zpos_;
    G4double holder_r_;
    G4double active_zpos_, cathode_zpos_, gate_zpos_, el_gap_zpos_, anode_zpos_;
    G4double gate_grid_zpos_, anode_grid_zpos_;


    // Vertex generators
    SegmentPointSampler* active_end_gen_;
    CylinderPointSampler2020* active_gen_;
    CylinderPointSampler2020* buffer_gen_;
    CylinderPointSampler2020* teflon_gen_;
    CylinderPointSampler2020* xenon_gen_;
    CylinderPointSampler2020* el_gap_gen_;
    CylinderPointSampler2020* hdpe_gen_;
    CylinderPointSampler2020* ring_gen_;
    CylinderPointSampler2020* cathode_gen_;
    CylinderPointSampler2020* gate_gen_;
    CylinderPointSampler2020* anode_gen_;
    CylinderPointSampler2020* holder_gen_;

    // Geometry Navigator
    G4Navigator* geom_navigator_;

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    // Logical volume where the class is placed
    G4LogicalVolume* mother_logic_;
    G4VPhysicalVolume* mother_phys_;
    G4VPhysicalVolume* active_phys_;
    G4VPhysicalVolume* buffer_phys_;
    G4Material* gas_;
    G4double pressure_;
    G4double temperature_;

    // Pointers to materials definition
    G4Material* hdpe_;
    G4Material* pe500_;
    G4Material* tpb_;
    G4Material* teflon_;
    G4Material* copper_;
    G4Material* steel_;

    G4double el_gap_gen_disk_diam_;
    G4double el_gap_gen_disk_x_, el_gap_gen_disk_y_;
    G4double el_gap_gen_disk_zmin_, el_gap_gen_disk_zmax_;
  };


  inline void Next100FieldCage::SetELtoSapphireWDWdistance(G4double distance){
    gate_sapphire_wdw_dist_ = distance;
  }

} //end namespace nexus
#endif
