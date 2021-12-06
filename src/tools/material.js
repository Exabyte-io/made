import math from "mathjs";
import {Basis} from "../basis/basis";
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

function scaleNonPeriodicLattice(material) {
    const lattice = new Lattice(material.lattice);
    const basis = new Basis({
        ...material.basis,
        cell: lattice.vectorArrays,
        units: "crystal"
    })
    basis.toCartesian();
    basis.toJSON();
    //TODO: Change to match the new scaling factor variable from SOF-5214
    const scalingFactor = 2.0;
    const newLattice = new Lattice({
        a: basis.maxPairwiseDistance * scalingFactor,
        type: "CUB"
    });
    newLattice.toJSON();
    return newLattice;
}

function getBasisConfigTranslatedToCenter(material) {
    const lattice = new Lattice(material.lattice);
    const basis = new Basis({
        ...material.basis,
        cell: lattice.vectorArrays
    });
    basis.toCartesian();
    const centerOfCoordinates = basis.centerOfCoordinatesPoint;
    const centerOfLattice = math.multiply(0.5, lattice.vectorArrays.reduce((a, b) => math.add(a, b)));
    const translationVector = math.subtract(centerOfLattice, centerOfCoordinates);
    basis.translateByVector(translationVector);
    const newBasisConfig = {
        elements: basis.elements,
        coordinates: basis.coordinates,
        cell: basis.cell,
        units: basis.units
    };
    return newBasisConfig;
}

export default {
    scaleOneLatticeVector,
    scaleNonPeriodicLattice,
    getBasisConfigTranslatedToCenter,
}
