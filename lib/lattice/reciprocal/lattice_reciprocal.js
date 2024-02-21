"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.ReciprocalLattice = void 0;
const constants_1 = require("@exabyte-io/code.js/dist/constants");
const array_almost_equal_1 = __importDefault(require("array-almost-equal"));
const lodash_1 = __importDefault(require("lodash"));
const math_1 = __importDefault(require("../../math"));
const lattice_1 = require("../lattice");
const paths_1 = require("./paths");
const symmetry_points_1 = require("./symmetry_points");
class ReciprocalLattice extends lattice_1.Lattice {
    /**
     * Get reciprocal vectors for the current Lattice in cartesian (2pi / a) units
     * @return {Array[]}
     */
    get reciprocalVectors() {
        const vectors_ = this.vectors.vectorArrays;
        const a = math_1.default.vlen(vectors_[0]);
        const divider = math_1.default.multiply(vectors_[0], math_1.default.cross(vectors_[1], vectors_[2])) / a;
        return [
            math_1.default.multiply(math_1.default.cross(vectors_[1], vectors_[2]), 1 / divider),
            math_1.default.multiply(math_1.default.cross(vectors_[2], vectors_[0]), 1 / divider),
            math_1.default.multiply(math_1.default.cross(vectors_[0], vectors_[1]), 1 / divider),
        ];
    }
    /**
     * Norms of reciprocal vectors.
     * @return {number[]}
     */
    get reciprocalVectorNorms() {
        return this.reciprocalVectors.map((vec) => math_1.default.norm(vec));
    }
    /**
     * Ratio of reciprocal vector norms scaled by the inverse of the largest component.
     * @return {number[]}
     */
    get reciprocalVectorRatios() {
        const norms = this.reciprocalVectorNorms;
        const maxNorm = math_1.default.max(...norms);
        return norms.map((n) => n / maxNorm);
    }
    /**
     * Get point (in crystal coordinates) in cartesian coordinates.
     * @param {Array} point - point in 3D space
     * @return {Array}
     */
    getCartesianCoordinates(point) {
        return math_1.default.multiply(point, this.reciprocalVectors);
    }
    /**
     * Get the list of high-symmetry points for the current lattice.
     * @return {Object[]}
     */
    get symmetryPoints() {
        return (0, symmetry_points_1.symmetryPoints)(this);
    }
    /**
     * Get the default path in reciprocal space for the current lattice.
     * @return {Array[]}
     */
    get defaultKpointPath() {
        return paths_1.paths[this.typeExtended] || paths_1.paths[this.type];
    }
    /**
     * Find/mark the high symmetry points on a list with raw data and return the edited list.
     * @param {Array} dataPoints - list of point coordinates
     * @return {Object[]}
     */
    extractKpointPath(dataPoints = []) {
        const kpointPath = [];
        const symmPoints = this.symmetryPoints;
        dataPoints.forEach((point, index) => {
            const symmPoint = symmPoints.find((x) => {
                return (0, array_almost_equal_1.default)(x.coordinates, point, 1e-4);
            });
            if (symmPoint) {
                kpointPath.push({
                    point: symmPoint.point,
                    steps: index,
                    coordinates: symmPoint.coordinates,
                });
            }
        });
        return kpointPath;
    }
    /**
     * Calculate grid dimension based on reciprocal lattice vectors.
     * @param {number} nPoints - Total number of points
     * @param {number} index - Index of reciprocal vector
     * @return {number} - Grid dimension in direction of reciprocal vector
     * @todo This could be moved to a separate KGrid class.
     */
    calculateDimension(nPoints, index) {
        const norms = this.reciprocalVectorNorms;
        const [j, k] = [0, 1, 2].filter((i) => i !== index); // get indices of other two dimensions
        const N = math_1.default.cbrt((nPoints * norms[index] ** 2) / (norms[j] * norms[k]));
        return math_1.default.max(1, math_1.default.ceil(N));
    }
    /**
     * Calculate grid dimensions from total number of k-points.
     * @param {number} nKpoints - Total number of k-points.
     * @return {number[]} - Grid dimensions
     */
    getDimensionsFromPointsCount(nKpoints) {
        const indices = [0, 1, 2];
        return indices.map((i) => this.calculateDimension(nKpoints, i));
    }
    get conversionTable() {
        const { a } = this;
        return {
            [constants_1.ATOMIC_COORD_UNITS.cartesian]: {
                [constants_1.units.angstrom]: (2 * math_1.default.PI) / a,
            },
            [constants_1.units.angstrom]: {
                [constants_1.ATOMIC_COORD_UNITS.cartesian]: a / (2 * math_1.default.PI),
            },
        };
    }
    /**
     * Calculate grid dimensions from k-point spacing, i.e.
     * the maximum distance between adjacent points along a reciprocal axis.
     * Note: just as the lattice vectors spacing is in cartesian (2pi / a) units by default
     * @param {number} spacing - maximum Spacing between k-points
     * @param {string} units - units of spacing parameter (default: 2pi / a)
     * @return {number[]}
     */
    getDimensionsFromSpacing(spacing, units = constants_1.ATOMIC_COORD_UNITS.cartesian) {
        const factor = this.conversionTable[units][constants_1.ATOMIC_COORD_UNITS.cartesian] || 1;
        return this.reciprocalVectorNorms.map((norm) => {
            return math_1.default.max(1, math_1.default.ceil(lodash_1.default.round(norm / (spacing * factor), 4)));
        });
    }
    /**
     * Calculate grid spacing as average of spacing along individual reciprocal axes.
     * @param {number[]} dimensions - Array of dimensions
     * @param {string} units - units of spacing parameter (default: 2pi / a)
     * @return {number} - average grid spacing
     */
    getSpacingFromDimensions(dimensions, units = constants_1.ATOMIC_COORD_UNITS.cartesian) {
        const factor = this.conversionTable[constants_1.ATOMIC_COORD_UNITS.cartesian][units] || 1;
        const norms = this.reciprocalVectorNorms;
        return factor * math_1.default.mean(dimensions.map((dim, i) => norms[i] / math_1.default.max(1, dim)));
    }
}
exports.ReciprocalLattice = ReciprocalLattice;
