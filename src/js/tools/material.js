import { math } from "@mat3ra/code/dist/js/math";

import { ATOMIC_COORD_UNITS } from "../constants";
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
    const originalUnits = material.Basis.units;
    material.toCartesian();
    const updatedBasis = material.Basis;
    const centerOfCoordinates = updatedBasis.centerOfCoordinatesPoint;
    const centerOfLattice = math.multiply(
        0.5,
        material.Lattice.vectorArrays.reduce((a, b) => math.add(a, b)),
    );
    const translationVector = math.subtract(centerOfLattice, centerOfCoordinates);
    updatedBasis.translateByVector(translationVector);
    material.setBasis(updatedBasis.toJSON());
    if (originalUnits !== ATOMIC_COORD_UNITS.cartesian) material.toCrystal();
}

export default {
    scaleOneLatticeVector,
    scaleLatticeToMakeNonPeriodic,
    getBasisConfigTranslatedToCenter,
};
