/*
 * @summary Scales one lattice vector for the given material
 * @param material {Material} The material acted upon.
 * @param key {String} Lattice vector key.
 * @param factor {Number} Float scaling factor.
 * @return {null}
 */

function scaleOneLatticeVector(material, key = 'a', factor = 1.0) {

    material.toCartesian();

    const lattice = material.lattice;
    lattice[key] = lattice[key] * factor;

    material.lattice = lattice;

    material.toCrystal();

}

export default {
    scaleOneLatticeVector,
}
