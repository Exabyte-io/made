declare namespace _default {
    export { latticePointsInSupercell };
}
export default _default;
/**
 * Returns the list of points on the original lattice contained in the supercell in fractional coordinates.
 * Source: https://pymatgen.org/_modules/pymatgen/util/coord.html
 */
declare function latticePointsInSupercell(supercellMatrix: any): import("../basis/types").Coordinate[];
