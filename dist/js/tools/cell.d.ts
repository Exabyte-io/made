import { Coordinate3DSchema, Matrix3X3Schema } from "@mat3ra/esse/dist/js/types";
/**
 * Returns the list of points on the original lattice contained in the supercell in fractional coordinates.
 * Source: https://pymatgen.org/_modules/pymatgen/util/coord.html
 */
declare function latticePointsInSupercell(supercellMatrix: Matrix3X3Schema): Coordinate3DSchema[];
declare const _default: {
    latticePointsInSupercell: typeof latticePointsInSupercell;
};
export default _default;
