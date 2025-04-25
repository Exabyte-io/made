import { LatticeSchema, Matrix3X3Schema } from "@mat3ra/esse/dist/js/types";
/**
 * Returns lattice vectors for a primitive cell for a lattice.
 * @param latticeConfig - Lattice config.
 * @return Cell.vectorsAsArray
 */
export declare function getPrimitiveLatticeVectorsFromConfig(latticeConfig: LatticeSchema): Matrix3X3Schema;
