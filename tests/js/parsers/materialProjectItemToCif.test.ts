import type { MaterialsProjectSchema } from "@mat3ra/esse/dist/js/types";
import { expect } from "chai";

import materialProjectItemToCif from "../../../src/js/parsers/materialProjectItemToCif";
import niobiumData from "../../fixtures/materials_project/mp-1094120.json";
import rubidiumData from "../../fixtures/materials_project/mp-1179802.json";
import carbonData from "../../fixtures/materials_project/mp-1197903.json";

type Schema = Pick<
    MaterialsProjectSchema,
    "material_id" | "formula_pretty" | "structure" | "symmetry" | "composition"
>;

describe("materialProjectItemToCif", () => {
    describe("Complete CIF Generation", () => {
        it("should generate complete CIF for carbon structure", () => {
            const cif = materialProjectItemToCif(carbonData as unknown as Schema);

            const expectedCif = `#\\#CIF_1.1
##########################################################################
#               Crystallographic Information Format file
#               Generated from Materials Project API
#
#  Material ID: mp-1197903
#  Formula: C
#  Generated from structure data
##########################################################################

data_mp_1197903
_symmetry_space_group_name_H-M          'P1'
_symmetry_Int_Tables_number             1
_symmetry_cell_setting                  Triclinic
_cell_length_a                          7.571152
_cell_length_b                          9.190471
_cell_length_c                          23.546004
_cell_angle_alpha                       97.892904
_cell_angle_beta                        96.750659
_cell_angle_gamma                       104.975528
_cell_volume                            1547.569425
_cell_formula_units_Z                   80
_chemical_formula_sum                  'C'
loop_
  _symmetry_equiv_pos_site_id
  _symmetry_equiv_pos_as_xyz
   1  'x, y, z'
loop_
  _atom_site_type_symbol
  _atom_site_label
  _atom_site_symmetry_multiplicity
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_attached_hydrogens
  _atom_site_B_iso_or_equiv
  _atom_site_occupancy
   C  C    1   0.274381   0.493366   0.004922 0  .  1.0
   C  C2   1   0.433610   0.480803   0.997758 0  .  1.0
   C  C3   1   0.577207   0.485377   0.036976 0  .  1.0
   C  C4   1   0.667848   0.519697   0.087519 0  .  1.0
   C  C5   1   0.675754   0.584601   0.145656 0  .  1.0
`;

            expect(cif).to.equal(expectedCif);
        });

        it("should generate complete CIF for niobium structure", () => {
            const cif = materialProjectItemToCif(niobiumData as unknown as Schema);

            const expectedCif = `#\\#CIF_1.1
##########################################################################
#               Crystallographic Information Format file
#               Generated from Materials Project API
#
#  Material ID: mp-1094120
#  Formula: Nb
#  Generated from structure data
##########################################################################

data_mp_1094120
_symmetry_space_group_name_H-M          'P1'
_symmetry_Int_Tables_number             1
_symmetry_cell_setting                  Triclinic
_cell_length_a                          5.374343
_cell_length_b                          6.938011
_cell_length_c                          7.179916
_cell_angle_alpha                       107.496791
_cell_angle_beta                        108.620266
_cell_angle_gamma                       102.624921
_cell_volume                            226.839615
_cell_formula_units_Z                   12
_chemical_formula_sum                  'Nb'
loop_
  _symmetry_equiv_pos_site_id
  _symmetry_equiv_pos_as_xyz
   1  'x, y, z'
loop_
  _atom_site_type_symbol
  _atom_site_label
  _atom_site_symmetry_multiplicity
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_attached_hydrogens
  _atom_site_B_iso_or_equiv
  _atom_site_occupancy
   Nb Nb   1   0.069490   0.140393   0.503316 0  .  1.0
   Nb Nb2  1   0.509667   0.212573   0.340606 0  .  1.0
   Nb Nb3  1   0.072719   0.731441   0.133756 0  .  1.0
   Nb Nb4  1   0.046112   0.405816   0.262338 0  .  1.0
   Nb Nb5  1   0.024193   0.983918   0.882255 0  .  1.0
   Nb Nb6  1   0.596368   0.925967   0.050041 0  .  1.0
   Nb Nb7  1   0.683872   0.720414   0.393695 0  .  1.0
   Nb Nb8  1   0.054694   0.564241   0.697763 0  .  1.0
   Nb Nb9  1   0.553340   0.563811   0.699148 0  .  1.0
   Nb Nb10 1   0.428080   0.395001   0.998307 0  .  1.0
   Nb Nb11 1   0.771339   0.251176   0.789395 0  .  1.0
   Nb Nb12 1   0.336504   0.880067   0.604147 0  .  1.0
`;

            expect(cif).to.equal(expectedCif);
        });

        it("should generate complete CIF for rubidium structure", () => {
            const cif = materialProjectItemToCif(rubidiumData as unknown as Schema);

            const expectedCif = `#\\#CIF_1.1
##########################################################################
#               Crystallographic Information Format file
#               Generated from Materials Project API
#
#  Material ID: mp-1179802
#  Formula: Rb
#  Generated from structure data
##########################################################################

data_mp_1179802
_symmetry_space_group_name_H-M          'P1'
_symmetry_Int_Tables_number             1
_symmetry_cell_setting                  Triclinic
_cell_length_a                          8.978624
_cell_length_b                          9.419573
_cell_length_c                          9.985111
_cell_angle_alpha                       102.415193
_cell_angle_beta                        114.107995
_cell_angle_gamma                       92.537353
_cell_volume                            744.494598
_cell_formula_units_Z                   8
_chemical_formula_sum                  'Rb'
loop_
  _symmetry_equiv_pos_site_id
  _symmetry_equiv_pos_as_xyz
   1  'x, y, z'
loop_
  _atom_site_type_symbol
  _atom_site_label
  _atom_site_symmetry_multiplicity
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_attached_hydrogens
  _atom_site_B_iso_or_equiv
  _atom_site_occupancy
   Rb Rb   1   0.388715   0.794630   0.234057 0  .  1.0
   Rb Rb2  1   0.974850   0.820843   0.782301 0  .  1.0
   Rb Rb3  1   0.894549   0.599739   0.269305 0  .  1.0
   Rb Rb4  1   0.780495   0.299401   0.750595 0  .  1.0
   Rb Rb5  1   0.273446   0.288511   0.733440 0  .  1.0
   Rb Rb6  1   0.453472   0.782270   0.726238 0  .  1.0
   Rb Rb7  1   0.849128   0.088731   0.229780 0  .  1.0
   Rb Rb8  1   0.364381   0.302596   0.249545 0  .  1.0
`;

            expect(cif).to.equal(expectedCif);
        });
    });
});
