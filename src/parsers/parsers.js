import cif from "./cif";
import espresso from "./espresso";
import poscar from "./poscar";
import { getLineFromContent } from "./utils";
import xyz from "./xyz";

/**
 * Function returns the of the file that should contain file contents based on the fileExtension value
 * @param {String} fileContent
 * @param {String} fileExtension
 * @param {String} failoverName
 * @returns {String}
 */
export function getNameFromContents(fileContent, fileExtension, failoverName = "material") {
    let lineNumber;
    if (fileExtension === "poscar") {
        lineNumber = 0;
    }
    if (fileExtension === "xyz") {
        lineNumber = 1;
    }
    const nameFromContent = getLineFromContent(fileContent, lineNumber);
    if (nameFromContent) return nameFromContent;
    return failoverName;
}

export default {
    xyz,
    poscar,
    cif,
    espresso,
    getNameFromContents,
};
