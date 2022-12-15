// ----------------------------------------------------------------------------
// nexus | SiPM11_eff.cc
//
// Geometry of a Hamamatsu 1x1 mm2 SiPM with 6x6 efficiencies
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "SiPM11_eff.h"
#include "SensorSD.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"

#include <G4Box.hh>
#include <G4GenericMessenger.hh>
#include <G4LogicalVolume.hh>
#include <G4VisAttributes.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4SDManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>

#include <CLHEP/Units/SystemOfUnits.h>


namespace nexus {

  using namespace CLHEP;

  SiPM11_eff::SiPM11_eff():
        GeometryBase(),
		    sipm_visibility_(true)

  {
    /// Messenger
    msg_ = new G4GenericMessenger(this, "/Geometry/Fib_SiPM_cyl/",
          "Control commands of geometry Fiber cylinder with SiPMs.");
    msg_->DeclareProperty("sipm_visibility", sipm_visibility_, "SiPM11_eff Visibility");
  }



  SiPM11_eff::~SiPM11_eff()
  {
    delete msg_;
  }



  G4ThreeVector SiPM11_eff::GetDimensions() const
  {
    return dimensions_;
  }



  void SiPM11_eff::Construct()
  {
    // PACKAGE ///////////////////////////////////////////////////////

    // G4double sipm_x = 2.425 * mm;
    // G4double sipm_y = 1.900 * mm;
    // G4double sipm_z = 0.850 * mm;
    G4double sipm_x = 1. * mm;
    G4double sipm_y = 1. * mm;
    G4double sipm_z = 0.850 * mm;

    dimensions_.setX(sipm_x);
    dimensions_.setY(sipm_y);
    dimensions_.setZ(sipm_z);

    G4Box* sipm_solid = new G4Box("SiPM11_eff", sipm_x/2., sipm_y/2., sipm_z/2);

    G4Material* epoxy = materials::Epoxy();
    epoxy->SetMaterialPropertiesTable(opticalprops::GlassEpoxy());

    G4LogicalVolume* sipm_logic =
      new G4LogicalVolume(sipm_solid, epoxy, "SiPM11_eff");

    this->SetLogicalVolume(sipm_logic);


    // // PCB ///////////////////////////////////////////////////////
    //
    // G4double pcb_z = 0.550 * mm;
    //
    // G4Material* plastic = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYCARBONATE");
    //
    // G4Box* plastic_solid = new G4Box("PLASTIC", sipm_x/2., sipm_y/2., pcb_z/2);
    //
    // G4LogicalVolume* plastic_logic =
    // new G4LogicalVolume(plastic_solid, plastic, "PLASTIC");
    //
    // G4double epoxy_z = 0.300 * mm;
    G4double epoxy_z = 0. * mm;
    //
    // new G4PVPlacement(0, G4ThreeVector(0, 0., epoxy_z/2), plastic_logic,
		//       "PLASTIC", sipm_logic, false, 0, false);

    // ACTIVE WINDOW /////////////////////////////////////////////////

    G4double active_side     = 1.0   * mm;
    G4double active_depth    = 0.01   * mm;
    // G4double active_offset_x = 0.975 * mm;
    G4double active_offset_x = 0. * mm;

    G4Box* active_solid =
      new G4Box("PHOTODIODES", active_side/2., active_side/2., active_depth/2);

    G4Material* silicon =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

    G4LogicalVolume* active_logic =
      new G4LogicalVolume(active_solid, silicon, "PHOTODIODES");

    G4double pos_x = - sipm_x/2. + active_offset_x + active_side/2.;
    // G4double pos_z = epoxy_z/2. - active_depth/2. - pcb_z/2;
    G4double pos_z = epoxy_z/2. + active_depth;

    new G4PVPlacement(0, G4ThreeVector(pos_x, 0., pos_z), active_logic,
		      "PHOTODIODES", sipm_logic, false, 0, false);


    // OPTICAL SURFACES //////////////////////////////////////////////

    // SiPM efficiency set using the official Hamamatsu specs.

    // const G4int entries = 26;
    const G4int entries = 22;

    // G4double energies[entries]     = {1.405*eV, 1.456*eV, 1.515*eV, 1.597*eV, 1.689*eV,
		// 		      1.763*eV, 1.836*eV, 1.915*eV, 2.007*eV, 2.105*eV,
		// 		      2.190*eV, 2.285*eV, 2.366*eV, 2.448*eV, 2.563*eV,
		// 		      2.718*eV, 2.838*eV, 2.977*eV, 3.099*eV, 3.243*eV,
		// 		      3.387*eV, 3.525*eV, 3.608*eV, 3.695*eV, 3.762*eV,
		// 		      3.857*eV };
    G4double energies[entries]     = {
      h_Planck * c_light / (864.707 * nm), h_Planck * c_light / (812.892 * nm),
      h_Planck * c_light / (766.143 * nm), h_Planck * c_light / (724.467 * nm),
      h_Planck * c_light / (686.591 * nm), h_Planck * c_light / (652.518 * nm),
      h_Planck * c_light / (627.310 * nm), h_Planck * c_light / (599.570 * nm),
      h_Planck * c_light / (573.096 * nm), h_Planck * c_light / (546.624 * nm),
      h_Planck * c_light / (520.150 * nm), h_Planck * c_light / (462.028 * nm),
      h_Planck * c_light / (423.991 * nm), h_Planck * c_light / (396.035 * nm),
      h_Planck * c_light / (381.996 * nm), h_Planck * c_light / (370.489 * nm),
      h_Planck * c_light / (357.717 * nm), h_Planck * c_light / (348.742 * nm),
      h_Planck * c_light / (343.562 * nm), h_Planck * c_light / (338.382 * nm),
      h_Planck * c_light / (331.937 * nm), h_Planck * c_light / (326.757 * nm)
    };
    G4double reflectivity[entries] = {0.      ,0.      ,0.      ,0.      ,0.      ,
				      0.      ,0.      ,0.      ,0.      ,0.      ,
				      0.      ,0.      ,0.      ,0.      ,0.      ,
              0.      ,0.      ,0.      ,0.      ,0.      ,
				      // 0.      ,0.      ,0.      ,0.      ,0.      ,
              // 0.      };
				      0.      ,0.       };
    G4double efficiency[entries]   = {0.0593, 0.0925, 0.1266, 0.1645, 0.2033,
      0.2445, 0.2873, 0.3300, 0.3728, 0.4163,
      0.4591, 0.5011, 0.4764, 0.4336, 0.3885,
      0.3433, 0.2981, 0.2529, 0.2070, 0.1610,
      0.1150, 0.0691
    };
    // G4double efficiency[entries]   = {0.0556  ,0.0698  ,0.0893  ,0.1250  ,0.1661  ,
		// 		      0.1983  ,0.2341  ,0.2663  ,0.3058  ,0.3488  ,
		// 		      0.3868  ,0.4247  ,0.4499  ,0.4734  ,0.4915  ,
		// 		      0.4991  ,0.4898  ,0.4662  ,0.4355  ,0.4002  ,
		// 		      0.3471  ,0.2878  ,0.2308  ,0.1620  ,0.0804  ,
		// 		      0.0390  };



    G4double efficiency_red[entries];
    for (G4int i=0; i<entries; ++i) {
      efficiency_red[i] = efficiency[i]*.6;
    }

    G4MaterialPropertiesTable* sipm_mt = new G4MaterialPropertiesTable();
    sipm_mt->AddProperty("EFFICIENCY", energies, efficiency_red, entries);
    sipm_mt->AddProperty("REFLECTIVITY", energies, reflectivity, entries);

    G4OpticalSurface* sipm_opsurf =
      new G4OpticalSurface("SIPM_OPSURF", unified, polished, dielectric_metal);
    sipm_opsurf->SetMaterialPropertiesTable(sipm_mt);

    new G4LogicalSkinSurface("SIPM_OPSURF", active_logic, sipm_opsurf);


    // SENSITIVE DETECTOR ////////////////////////////////////////////

    G4String sdname = "/SiPM11_eff/SiPM";
    G4SDManager* sdmgr = G4SDManager::GetSDMpointer();

    if (!sdmgr->FindSensitiveDetector(sdname, false)) {
      SensorSD* sipmsd = new SensorSD(sdname);
      sipmsd->SetDetectorVolumeDepth(1);
      sipmsd->SetDetectorNamingOrder(1000.);
      sipmsd->SetTimeBinning(1.*microsecond);
      sipmsd->SetMotherVolumeDepth(2);

      G4SDManager::GetSDMpointer()->AddNewDetector(sipmsd);
      active_logic->SetSensitiveDetector(sipmsd);
    }

    // Visibilities
    if (sipm_visibility_) {
      G4VisAttributes sipm_col = nexus::DirtyWhite();
      sipm_logic->SetVisAttributes(sipm_col);
      G4VisAttributes blue_col = nexus::Blue();
      blue_col.SetForceSolid(true);
      active_logic->SetVisAttributes(blue_col);
      G4VisAttributes plastic_col = nexus::Lilla();
      plastic_col.SetForceSolid(true);
      // plastic_logic->SetVisAttributes(plastic_col);
    }
    else {
      sipm_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      active_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      // plastic_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }
  }


} // end namespace nexus
