import type { MaterialsProjectSchema } from "@mat3ra/esse/dist/js/types";
type Schema = Pick<MaterialsProjectSchema, "material_id" | "formula_pretty" | "structure" | "symmetry" | "composition">;
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
export default function materialProjectItemToCif(materialProjectItem: Schema): string;
export {};
