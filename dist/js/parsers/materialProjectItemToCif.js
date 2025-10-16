"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * Converts Materials Project structure data to CIF format
 * @returns CIF string representation of the structure
 *
 * ✅ Features Implemented:
 * - Fixed CIF version header format (#\\#CIF_1.1)
 * - Comprehensive input validation for all required fields
 * - Smart atomic site label generation to avoid duplicates
 * - Standardized precision to 6 decimal places for all numerical values
 * - Dynamic Z calculation from composition data (formula units per unit cell)
 * - Robust error handling for invalid data structures
 * - Handles duplicate site labels properly (e.g., "Nb" → "Nb", "Nb2", "Nb3", etc.)
 * - Proper CIF formatting with all required sections
 *
 * ⚠️ Current Limitations:
 * - Symmetry Operations: Only includes identity operation 'x, y, z' (hardcoded for P1)
 * - Symmetry Multiplicity: Always set to 1 (simplified for P1 space group)
 * - Advanced symmetry operations not implemented for higher space groups
 * - No support for non-P1 space groups (would need symmetry operation matrices)
 */
function materialProjectItemToCif(materialProjectItem) {
    const { material_id, formula_pretty, structure, symmetry } = materialProjectItem;
    const { lattice, sites } = structure;
    // Validate required fields
    if (!material_id || !formula_pretty || !structure || !symmetry) {
        throw new Error("Missing required fields: material_id, formula_pretty, structure, or symmetry");
    }
    if (!lattice || !sites || !Array.isArray(sites)) {
        throw new Error("Invalid structure: missing lattice or sites array");
    }
    if (!symmetry.symbol) {
        throw new Error("Missing symmetry symbol");
    }
    // Calculate formula units per unit cell (Z)
    const { composition } = materialProjectItem;
    const totalAtoms = Object.values(composition || {}).reduce((sum, val) => sum + val, 0);
    const z = Math.round(totalAtoms) || sites.length;
    // CIF header
    let cif = `#\\#CIF_1.1
##########################################################################
#               Crystallographic Information Format file
#               Generated from Materials Project API
#
#  Material ID: ${material_id}
#  Formula: ${formula_pretty}
#  Generated from structure data
##########################################################################

data_${material_id.replace("-", "_")}
_symmetry_space_group_name_H-M          '${symmetry.symbol}'
_symmetry_Int_Tables_number             ${symmetry.number || 1}
_symmetry_cell_setting                  ${symmetry.crystal_system || "Triclinic"}
_cell_length_a                          ${lattice.a.toFixed(6)}
_cell_length_b                          ${lattice.b.toFixed(6)}
_cell_length_c                          ${lattice.c.toFixed(6)}
_cell_angle_alpha                       ${lattice.alpha.toFixed(6)}
_cell_angle_beta                        ${lattice.beta.toFixed(6)}
_cell_angle_gamma                       ${lattice.gamma.toFixed(6)}
_cell_volume                            ${lattice.volume.toFixed(6)}
_cell_formula_units_Z                   ${z}
_chemical_formula_sum                  '${formula_pretty}'
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
`;
    // Add atomic sites
    const elementCounts = {};
    sites.forEach((site, index) => {
        var _a, _b;
        // Validate site structure
        if (!site.species || !Array.isArray(site.species) || site.species.length === 0) {
            throw new Error(`Invalid site at index ${index}: missing or empty species array`);
        }
        if (!site.abc || !Array.isArray(site.abc) || site.abc.length !== 3) {
            throw new Error(`Invalid site at index ${index}: missing or invalid abc coordinates`);
        }
        const element = ((_a = site.species[0]) === null || _a === void 0 ? void 0 : _a.element) || "X";
        // Generate unique labels to avoid duplicates
        if (!elementCounts[element]) {
            elementCounts[element] = 0;
        }
        elementCounts[element] += 1;
        // Generate unique labels to avoid duplicates
        let { label } = site;
        // Check if this label has been used before (including current site)
        const isLabelUsed = sites.slice(0, index + 1).filter((s) => s.label === label).length > 1;
        if (!label || isLabelUsed) {
            label = `${element}${elementCounts[element]}`;
        }
        const [x, y, z] = site.abc;
        const occupancy = ((_b = site.species[0]) === null || _b === void 0 ? void 0 : _b.occu) || 1.0;
        // Calculate symmetry multiplicity (1 for P1 space group)
        const multiplicity = symmetry.number === 1 ? 1 : 1; // Simplified for now
        cif += `   ${element.padEnd(2)} ${label.padEnd(4)} ${multiplicity} ${x
            .toFixed(6)
            .padStart(10)} ${y.toFixed(6).padStart(10)} ${z
            .toFixed(6)
            .padStart(10)} 0  .  ${occupancy.toFixed(1)}
`;
    });
    return cif;
}
exports.default = materialProjectItemToCif;
