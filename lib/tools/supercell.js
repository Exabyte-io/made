"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const lattice_bravais_1 = require("../lattice/lattice_bravais");
const math_1 = __importDefault(require("../math"));
const cell_1 = __importDefault(require("./cell"));
const ADD = math_1.default.add;
/**
 * @summary Generates new basis for a supercell. For each site from basis generates shifts that are within supercell.
 */
function generateNewBasisWithinSupercell(basis, cell, supercell, supercellMatrix) {
    const oldBasis = basis.clone();
    const newBasis = basis.clone({ isEmpty: true });
    oldBasis.toCrystal();
    newBasis.toCrystal();
    oldBasis.elements.forEach((element) => {
        const coordinate = oldBasis.getCoordinateByIndex(element.id);
        const cartesianCoordinate = cell.convertPointToCartesian(coordinate);
        const shifts = cell_1.default.latticePointsInSupercell(supercellMatrix);
        shifts.forEach((comb) => {
            // "combination" is effectively a point in fractional coordinates here, hence the below
            const newPoint = ADD(cartesianCoordinate, supercell.convertPointToCartesian(comb));
            if (supercell.isPointInsideCell(newPoint)) {
                newBasis.addAtom({
                    element: element.value,
                    coordinate: supercell.convertPointToFractional(newPoint),
                });
            }
        });
    });
    return newBasis;
}
/**
 * @summary Generates supercell config for the specified material.
 * @param material
 * @param supercellMatrix {Number[][]}
 */
function generateConfig(material, supercellMatrix) {
    const det = math_1.default.det(supercellMatrix);
    if (det === 0) {
        throw new Error("Scaling matrix is degenerate.");
    }
    const cell = material.Lattice.Cell;
    const supercell = cell.cloneAndScaleByMatrix(supercellMatrix);
    const newBasis = generateNewBasisWithinSupercell(material.Basis, cell, supercell, supercellMatrix);
    const newLattice = lattice_bravais_1.LatticeBravais.fromVectors({
        a: supercell.vector1,
        b: supercell.vector2,
        c: supercell.vector3,
    });
    return {
        name: `${material.name} - supercell ${JSON.stringify(supercellMatrix)}`,
        basis: newBasis.toJSON(),
        lattice: newLattice.toJSON(),
    };
}
exports.default = {
    generateConfig,
    generateNewBasisWithinSupercell,
};
