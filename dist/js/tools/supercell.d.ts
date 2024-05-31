import { Basis } from "../basis/basis";
import { Cell } from "../cell/cell";
import { Material } from "../material";
/**
 * @summary Generates new basis for a supercell. For each site from basis generates shifts that are within supercell.
 */
declare function generateNewBasisWithinSupercell(basis: Basis, cell: Cell, supercell: Cell, supercellMatrix: number[][]): Basis;
/**
 * @summary Generates supercell config for the specified material.
 * @param material
 * @param supercellMatrix {Number[][]}
 */
declare function generateConfig(material: Material, supercellMatrix: number[][]): {
    name: string;
    basis: import("../basis/basis").BasisSchema;
    lattice: import("@mat3ra/esse/lib/js/types").LatticeImplicitSchema;
};
declare const _default: {
    generateConfig: typeof generateConfig;
    generateNewBasisWithinSupercell: typeof generateNewBasisWithinSupercell;
};
export default _default;
