import cif from "./cif";
// eslint-disable-next-line import/no-cycle
import espresso from "./espresso";
import poscar from "./poscar";
// eslint-disable-next-line import/no-cycle
import xyz from "./xyz";

/**
 * Function returns the number of atoms in a file using the proper parser function based on the file extension.
 * @param {String} fileExtension
 * @param {String} fileContent
 * @returns {Number}
 */
export function getNumberOfAtomsInFileByExtension(fileExtension, fileContent) {
    let numberOfAtoms = 0;
    if (fileExtension === "poscar") {
        numberOfAtoms = poscar.poscarFileAtomsCount(fileContent);
    }
    if (fileExtension === "xyz") {
        numberOfAtoms = xyz.xyzFileAtomsCount(fileContent);
    }
    return numberOfAtoms;
}

export default {
    xyz,
    poscar,
    cif,
    espresso,
    getNumberOfAtomsInFileByExtension,
};
