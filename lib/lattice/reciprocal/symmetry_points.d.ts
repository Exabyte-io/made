import { Lattice } from "../lattice";
/**
 * Returns a list of symmetry points for the specified lattice.
 */
export declare function symmetryPoints(lattice: Lattice): {
    point: string;
    coordinates: number[];
}[];
