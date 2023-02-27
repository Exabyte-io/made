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
     * @param nPoints - Total number of points
     * @param index - Index of reciprocal vector
     * @return {number} - Grid dimension in direction of reciprocal vector
     * @todo This could be moved to a separate KGrid class.
     */
    calculateDimension(nPoints, index) {
        const norms = this.reciprocalVectorNorms;
        const [j, k] = [0, 1, 2].filter((i) => i !== index); // get indices of other two dimensions
        const N = Math.cbrt((nPoints * norms[index] ** 2) / (norms[j] * norms[k]));
        return Math.max(1, Math.ceil(N));
    }

    /**
     * Calculate grid dimensions from total number of k-points.
     * @param nKpoints - Total number of k-points.
     * @return {number[]} - Grid dimensions
     */
    getDimensionsFromPoints(nKpoints) {
        const indices = [0, 1, 2];
        return indices.map((i) => this.calculateDimension(nKpoints, i));
    }

    /**
     * Calculate grid dimensions from k-point spacing, i.e.
     * the maximum distance between adjacent points along a reciprocal axis.
     * @param {number} spacing - maximum Spacing between k-points (1/Ang)
     * @return {number[]}
     */
    getDimensionsFromSpacing(spacing) {
        return this.reciprocalVectorNorms.map((norm) => {
            return Math.max(1, Math.ceil(lodash.round(norm / spacing, 4)));
        });
    }
}
