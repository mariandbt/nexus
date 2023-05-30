#include "Fiber_barrel_meth.h"
#include "GenericWLSFiber.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
#include "CylinderPointSampler.h"
#include "FactoryBase.h"

#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4NistManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>


using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(Fiber_barrel_meth,GeometryBase)

Fiber_barrel_meth::Fiber_barrel_meth():GeometryBase(), radius_(1.*mm), length_(1.*cm),  radius_cyl_(1. *cm), cyl_vertex_gen_(0){
    msg_=new G4GenericMessenger(this,"/Geometry/Fiber_barrel_meth/","Control commands of geometry OpticalFibre");

    G4GenericMessenger::Command& radius_cmd =
            msg_->DeclareProperty("radius",radius_,"Radius of the cylindrical optical fibre");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius",false);
    radius_cmd.SetRange("radius>0.");

    G4GenericMessenger::Command& length_cmd =
            msg_->DeclareProperty("length",length_,"Length of the cylindrical optical fibre");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length",false);
    length_cmd.SetRange("length>0.");

    G4GenericMessenger::Command& radius_cyl_cmd =
            msg_->DeclareProperty("radius_cyl",radius_cyl_,"Radius of the cylinder");
    radius_cyl_cmd.SetUnitCategory("Length");
    radius_cyl_cmd.SetParameterName("radius_cyl",false);
    radius_cyl_cmd.SetRange("radius_cyl>0.");

    cyl_vertex_gen_ = new CylinderPointSampler(0.5*radius_cyl_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);
}
Fiber_barrel_meth::~Fiber_barrel_meth() {
    delete cyl_vertex_gen_;
    delete msg_;
}
void Fiber_barrel_meth::Construct(){
    // G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);
    G4Box* lab_solid = new G4Box("LAB", 20 * mm,20 * mm,10.1*cm);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);


    G4double x;
    G4double y;
    G4double z = 0.;
    G4double rot_angle;

    // Outer teflon cilynder

    G4String teflon_name = "Teflon_panel";

    G4double teflon_thickness = .1 * mm;

    G4Material* teflon_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
    teflon_mat->SetMaterialPropertiesTable(opticalprops::PTFE());

    G4Tubs* teflon_solid_vol =
      new G4Tubs(teflon_name, radius_cyl_ + radius_, radius_cyl_ + radius_ + teflon_thickness, length_/2., 0., 360.*deg);

    G4LogicalVolume* teflon_logic_vol =
      new G4LogicalVolume(teflon_solid_vol, teflon_mat, teflon_name);

    G4OpticalSurface* teflon_opsur =
    new G4OpticalSurface("teflon_OPSURF", unified, polished, dielectric_metal);
    teflon_opsur->SetMaterialPropertiesTable(opticalprops::PTFE());
    new G4LogicalSkinSurface("teflon_OPSURF", teflon_logic_vol, teflon_opsur);

    G4VisAttributes teflon_col = nexus::Lilla();
    teflon_logic_vol->SetVisAttributes(teflon_col);
    // teflon_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4ThreeVector teflon_pos = G4ThreeVector(0., 0., 0.);
    G4RotationMatrix* teflon_rot = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    teflon_rot->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*teflon_rot, teflon_pos),
                      teflon_logic_vol, teflon_name, lab_logic,
                      false, 0, false);



    // Inner methacrylate cilynder
    G4Material* tpb = materials::TPB();

    G4String window_name = "Methacrylate window";

    G4double window_thickness = 1.5 * mm;

    G4Material* window_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
    window_mat->SetMaterialPropertiesTable(opticalprops::PMMA());

    G4Tubs* window_solid_vol =
    new G4Tubs(window_name, radius_cyl_ - radius_ - window_thickness, radius_cyl_ - radius_, length_/2., 0., 360.*deg);

    G4LogicalVolume* window_logic_vol =
      new G4LogicalVolume(window_solid_vol, window_mat, window_name);
    G4VisAttributes window_col = nexus::LightBlue();
    window_logic_vol->SetVisAttributes(window_col);

    G4OpticalSurface* window_opsur =
      new G4OpticalSurface("window_OPSURF", unified, polished, dielectric_dielectric);
      // window_opsur->SetMaterialPropertiesTable(opticalprops::PMMA());
      window_opsur->SetMaterialPropertiesTable(opticalprops::TPB());
    new G4LogicalSkinSurface("window_OPSURF", window_logic_vol, window_opsur);

    G4ThreeVector window_pos = G4ThreeVector(0., 0., 0.);

    G4RotationMatrix* window_rot_ = new G4RotationMatrix();
    // rot_angle = pi/2.;
    rot_angle = 0.;
    window_rot_->rotateY(rot_angle);
    new G4PVPlacement(G4Transform3D(*window_rot_, window_pos),
                      window_logic_vol, window_name, lab_logic,
                      false, 0, false);


    // Fibers

    G4double n_fibers = radius_cyl_*pi/radius_;
    G4double dif_theta = .5*2*pi/n_fibers; // angular separation between fibers
    G4double theta;

    G4Material* ps = materials::PS();

    GenericWLSFiber* fiber_;
    G4LogicalVolume* fiber_logic;

    for (int i=0; i<=2*n_fibers; i++){

      fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
      fiber_->SetCoreOpticalProperties(opticalprops::Y11());
      fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
      fiber_->Construct();
      fiber_logic = fiber_->GetLogicalVolume();

      theta = dif_theta*i;
      x = radius_cyl_ * cos(theta);
      y = radius_cyl_ * sin(theta);
      new G4PVPlacement(0,G4ThreeVector(x, y, z),fiber_logic,
                              fiber_logic->GetName(),lab_logic,true,0,true);
    }


}


G4ThreeVector Fiber_barrel_meth::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);

    // // G4ThreeVector vertex(1.,1.,1.);
    // G4ThreeVector vertex(0., 0., 0.);
    // return vertex;

}
