"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
const cell_1 = require("../cell/cell");
const math_1 = __importDefault(require("../math"));
/**
 * Returns the list of points on the original lattice contained in the supercell in fractional coordinates.
 * Source: https://pymatgen.org/_modules/pymatgen/util/coord.html
 */
function latticePointsInSupercell(supercellMatrix) {
    const supercell = new cell_1.Cell(supercellMatrix);
    const diagonals = [
        [0, 0, 0],
        [0, 0, 1],
        [0, 1, 0],
        [0, 1, 1],
        [1, 0, 0],
        [1, 0, 1],
        [1, 1, 0],
        [1, 1, 1],
    ];
    const d_points = diagonals.map((p) => supercell.convertPointToCartesian(p));
    const mins = [0, 1, 2].map((i) => math_1.default.min(...d_points.map((p) => p[i])));
    const maxes = [0, 1, 2].map((i) => math_1.default.max(...d_points.map((p) => p[i])) + 1);
    const points = [];
    for (let i = mins[0]; i <= maxes[0]; i++) {
        for (let j = mins[1]; j <= maxes[1]; j++) {
            for (let k = mins[2]; k <= maxes[2]; k++) {
                points.push(supercell.convertPointToFractional([i, j, k]));
            }
        }
    }
    return points;
}
exports.default = {
    latticePointsInSupercell,
};
