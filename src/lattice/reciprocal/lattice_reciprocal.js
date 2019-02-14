import almostEqual from "array-almost-equal";
import math from "../../math";
import {Lattice} from "../lattice";
import {paths} from "./paths";
import {symmetryPoints} from "./symmetry_points";

export class ReciprocalLattice extends Lattice {

    /**
     * Create a Reciprocal lattice class.
     * @param {Object} config - same as for Lattice.
     */
    constructor(config) {
        super(config);
    }

    /**
     * Ger reciprocal vectors for the current Lattice.
     * @return {Array[]}
     */
    get reciprocalVectors() {
        const vectors_ = this.vectors.vectorArrays;
        const divider = math.multiply(vectors_[0], math.cross(vectors_[1], vectors_[2]));
        return [
            math.multiply(math.cross(vectors_[1], vectors_[2]), 1 / divider),
            math.multiply(math.cross(vectors_[0], vectors_[2]), 1 / divider),
            math.multiply(math.cross(vectors_[1], vectors_[0]), 1 / divider)
        ];
    }

    /**
     * Get point (in crystal coordinates) in cartesian coordinates.
     * @param {Array} point - point in 3D space
     * @return {Array}
     */
    getCartesianCoordinates(point) {return math.multiply(point, math.inv(this.reciprocalVectors))}

    /**
     * Get the list of high-symmetry points for the current lattice.
     * @return {Object[]}
     */
    get symmetryPoints() {return symmetryPoints(this)}

    /**
     * Get the default path in reciprocal space for the current lattice.
     * @return {Array[]}
     */
    get defaultKpointPath() {return paths[this.typeExtended] || paths[this.type]}

    /**
     * Find/mark the high symmetry points on a list with raw data and return the edited list.
     * @param {Array} dataPoints - list of point coordinates
     * @return {Object[]}
     */
    extractKpointPath(dataPoints = []) {
        const kpointPath = [];
        const symmPoints = this.symmetryPoints;

        dataPoints.forEach((point, index) => {
            const symmPoint = symmPoints.find(x => {
                return almostEqual(x.coordinates, point, 1E-4);
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

}
