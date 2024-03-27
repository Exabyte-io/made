import { ATOMIC_COORD_UNITS, units as UNITS } from "@mat3ra/code/dist/js/constants";
import almostEqual from "array-almost-equal";
import lodash from "lodash";

import math from "../../math";
import { Lattice } from "../lattice";
import { paths } from "./paths";
import { symmetryPoints } from "./symmetry_points";

export class ReciprocalLattice extends Lattice {
    /**
     * Get reciprocal vectors for the current Lattice in cartesian (2pi / a) units
     * @return {Array[]}
     */
    get reciprocalVectors() {
        const vectors_ = this.vectors.vectorArrays;
        const a = math.vlen(vectors_[0]);
        const divider = math.multiply(vectors_[0], math.cross(vectors_[1], vectors_[2])) / a;
        return [
            math.multiply(math.cross(vectors_[1], vectors_[2]), 1 / divider),
            math.multiply(math.cross(vectors_[2], vectors_[0]), 1 / divider),
            math.multiply(math.cross(vectors_[0], vectors_[1]), 1 / divider),
        ];
    }

    /**
     * Norms of reciprocal vectors.
     * @return {number[]}
     */
    get reciprocalVectorNorms() {
        return this.reciprocalVectors.map((vec) => math.norm(vec));
    }

    /**
     * Ratio of reciprocal vector norms scaled by the inverse of the largest component.
     * @return {number[]}
     */
    get reciprocalVectorRatios() {
        const norms = this.reciprocalVectorNorms;
        const maxNorm = math.max(...norms);
        return norms.map((n) => n / maxNorm);
    }

    /**
     * Get point (in crystal coordinates) in cartesian coordinates.
     * @param {Array} point - point in 3D space
     * @return {Array}
     */
    getCartesianCoordinates(point) {
        return math.multiply(point, this.reciprocalVectors);
    }

    /**
     * Get the list of high-symmetry points for the current lattice.
     * @return {Object[]}
     */
    get symmetryPoints() {
        return symmetryPoints(this);
    }

    /**
     * Get the default path in reciprocal space for the current lattice.
     * @return {Array[]}
     */
    get defaultKpointPath() {
        return paths[this.typeExtended] || paths[this.type];
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
                return almostEqual(x.coordinates, point, 1e-4);
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
        const N = math.cbrt((nPoints * norms[index] ** 2) / (norms[j] * norms[k]));
        return math.max(1, math.ceil(N));
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
            [ATOMIC_COORD_UNITS.cartesian]: {
                [UNITS.angstrom]: (2 * math.PI) / a,
            },
            [UNITS.angstrom]: {
                [ATOMIC_COORD_UNITS.cartesian]: a / (2 * math.PI),
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
    getDimensionsFromSpacing(spacing, units = ATOMIC_COORD_UNITS.cartesian) {
        const factor = this.conversionTable[units][ATOMIC_COORD_UNITS.cartesian] || 1;
        return this.reciprocalVectorNorms.map((norm) => {
            return math.max(1, math.ceil(lodash.round(norm / (spacing * factor), 4)));
        });
    }

    /**
     * Calculate grid spacing as average of spacing along individual reciprocal axes.
     * @param {number[]} dimensions - Array of dimensions
     * @param {string} units - units of spacing parameter (default: 2pi / a)
     * @return {number} - average grid spacing
     */
    getSpacingFromDimensions(dimensions, units = ATOMIC_COORD_UNITS.cartesian) {
        const factor = this.conversionTable[ATOMIC_COORD_UNITS.cartesian][units] || 1;
        const norms = this.reciprocalVectorNorms;
        return factor * math.mean(dimensions.map((dim, i) => norms[i] / math.max(1, dim)));
    }
}
