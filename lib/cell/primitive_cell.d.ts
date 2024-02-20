import { LatticeImplicitSchema } from "@mat3ra/esse/lib/js/types";
import { VectorsAsArray } from "../lattice/types";
/**
 * Returns lattice vectors for a primitive cell for a lattice.
 * @param lattice - Lattice instance.
 * @param  skipRounding - whether to skip rounding the lattice vectors.
 * @return Cell.vectorsAsArray
 */
export declare function primitiveCell(lattice: LatticeImplicitSchema, skipRounding?: boolean): VectorsAsArray;
