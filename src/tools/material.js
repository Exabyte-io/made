import {Basis} from "../basis/basis";
import {Lattice, nonPeriodicLatticeScalingFactor} from "../lattice/lattice";

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
 * Returns the lattice for a non-periodic structure. The lattice has been scaled to accomodate the size of the
 * non-periodic structure.
 * @param material {Material}
 * @returns {Object}
 */
function scaleLatticeToMakeNonPeriodic(material) {
    const basis = new Basis({
        ...material.Basis.toJSON(),
    })
    basis.toCartesian();
    const scalingFactor = nonPeriodicLatticeScalingFactor;
    const newLattice = new Lattice({
        a: basis.maxPairwiseDistance * scalingFactor,
        type: "CUB"
    });

    return newLattice.toJSON();
}

export default {
    scaleOneLatticeVector,
    scaleLatticeToMakeNonPeriodic,
}
