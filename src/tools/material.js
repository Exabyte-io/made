import math from "mathjs";

import { Lattice } from "../lattice/lattice";

/**
 * Scales one lattice vector for the given material
 * @param material {Material} The material acted upon.
 * @param key {String} Lattice vector key.
 * @param factor {Number} Float scaling factor.
 */

function scaleOneLatticeVector(material, key = "a", factor = 1.0) {
    material.toCartesian();

    const { lattice } = material;
    lattice[key] *= factor;

    material.lattice = lattice;

    material.toCrystal();
}

/**
 * Updates the size of a materials lattice using the minimumLatticeSize function.
 * The new size of the material is calculated based on the materials basis.
 * @param material {Material}
 */
function scaleLatticeToMakeNonPeriodic(material) {
    material.toCartesian();
    material.lattice = new Lattice({
        a: material.Basis.getMinimumLatticeSize(),
        type: "CUB",
    });
}

/**
 * Updates the basis of a material by translating the coordinates
 * so that the center of the material and lattice are aligned.
 * @param material {Material}
 * */
function getBasisConfigTranslatedToCenter(material) {
    const updatedBasis = material.Basis;
    const centerOfCoordinates = updatedBasis.centerOfCoordinatesPoint;
    const centerOfLattice = math.multiply(
        0.5,
        material.lattice.vectorArrays.reduce((a, b) => math.add(a, b)),
    );
    const translationVector = math.subtract(centerOfLattice, centerOfCoordinates);
    updatedBasis.translateByVector(translationVector);
    material.setBasis(updatedBasis.toJSON());
}

export default {
    scaleOneLatticeVector,
    scaleLatticeToMakeNonPeriodic,
    getBasisConfigTranslatedToCenter,
};
