import {Lattice} from "../lattice/lattice";

/**
 * Scales one lattice vector for the given material
 * @param material {Material} The material acted upon.
 * @param key {String} Lattice vector key.
 * @param factor {Number} Float scaling factor.
 */

function scaleOneLatticeVector(material, key = 'a', factor = 1.0) {

    material.toCartesian();

    const lattice = material.lattice;
    lattice[key] = lattice[key] * factor;

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
        type: "CUB"
    });
}

export default {
    scaleOneLatticeVector,
    scaleLatticeToMakeNonPeriodic,
}
