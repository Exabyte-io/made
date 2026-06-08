import { MaterialSchema } from "@mat3ra/esse/dist/js/types";
/**
 * Strip VASP-style comments (everything after !) from POSCAR content.
 * @param poscarContent - POSCAR file content with potential comments.
 * @return POSCAR content without comments.
 */
declare function stripPoscarComments(poscarContent: string): string;
/**
 * Obtain a textual representation of a material in POSCAR format.
 * @param materialOrConfig - material class instance or config object.
 * @param omitConstraints - whether to discard constraints passed with material.
 */
declare function toPoscar(materialOrConfig: MaterialSchema, omitConstraints?: boolean): string;
/**
 * @summary calculates the number of atoms in a poscar file based on the summation of the numbers in line 7 of the file.
 * Poscar file formatting: https://www.vasp.at/wiki/index.php/POSCAR
 */
export declare function atomsCount(poscarFileContent: string): number;
/**
 * Parses POSCAR file into a Material config object.
 * @param fileContent - POSCAR file content.
 * @return Material config.
 */
declare function fromPoscar(fileContent: string): object;
/**
 * @summary Checks if a string has a POSCAR format (first 8 lines are read)
 * @param text - string to check
 */
declare function isPoscar(text: string): boolean;
declare const _default: {
    isPoscar: typeof isPoscar;
    toPoscar: typeof toPoscar;
    fromPoscar: typeof fromPoscar;
    atomicConstraintsCharFromBool: (bool: boolean) => string;
    atomsCount: typeof atomsCount;
    stripPoscarComments: typeof stripPoscarComments;
};
export default _default;
