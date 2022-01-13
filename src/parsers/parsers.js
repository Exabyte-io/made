import { defaultMaterialConfig } from "../default_material";
import cif from "./cif";
import espresso from "./espresso";
import poscar from "./poscar";
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
        numberOfAtoms = poscar.getAtomsCount(fileContent);
    }
    if (fileExtension === "xyz") {
        numberOfAtoms = xyz.getAtomsCount(fileContent);
    }
    return numberOfAtoms;
}

/**
 * Function converts an XYZ formatted structure as a POSCAR formatted structure
 * @param {String} xyzContent
 * @returns {String}
 */
export function xyzToPoscar(xyzContent) {
    const xyzConfig = defaultMaterialConfig;
    const xyzArray = xyzContent.split(/\r?\n/);
    const xyzArrayBasisOnly = xyzArray.slice(2, -1);
    const xyzBasis = xyzArrayBasisOnly.join("\n");
    xyzConfig.basis = xyz.toBasisConfig(xyzBasis);
    xyzConfig.basis.units = "cartesian";
    return poscar.toPoscar(xyzConfig);
}

export default {
    xyz,
    poscar,
    cif,
    espresso,
    getNumberOfAtomsInFileByExtension,
    xyzToPoscar,
};
