"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const math_1 = require("@mat3ra/code/dist/js/math");
const constants_1 = require("../constants");
const lattice_1 = require("../lattice/lattice");
/**
 * Scales one lattice vector for the given material
 * @param material {Material} The material acted upon.
 * @param key {String} Lattice vector key.
 * @param factor {Number} Float scaling factor.
 */
function scaleOneLatticeVector(material, key = "a", factor = 1.0) {
    var _a;
    material.toCartesian();
    const { lattice } = material;
    // @ts-ignore
    (_a = lattice.vectors) === null || _a === void 0 ? void 0 : _a[key] = lattice.vectors[key].map((v) => v * factor);
    material.lattice = lattice;
    material.toCrystal();
}
/**
 * Updates the size of a materials lattice using the minimumLatticeSize function.
 * The new size of the material is calculated based on the materials basis.
 * @param material {Material}
 */
function scaleLatticeToMakeNonPeriodic(material) {
    // @ts-ignore
    material.lattice = lattice_1.Lattice.fromConfigPartial({
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
    const centerOfLattice = math_1.math.multiply(0.5, material.Lattice.vectorArrays.reduce((a, b) => math_1.math.add(a, b)));
    const translationVector = math_1.math.subtract(centerOfLattice, centerOfCoordinates);
    updatedBasis.translateByVector(translationVector);
    material.setBasis(updatedBasis.toJSON());
    if (originalUnits !== constants_1.ATOMIC_COORD_UNITS.cartesian)
        material.toCrystal();
}
exports.default = {
    scaleOneLatticeVector,
    scaleLatticeToMakeNonPeriodic,
    getBasisConfigTranslatedToCenter,
};
