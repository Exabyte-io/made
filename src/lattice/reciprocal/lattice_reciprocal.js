import almostEqual from "array-almost-equal";
import math from "../../math";
import {Lattice} from "../lattice";
import {paths} from "./paths";
import {symmetryPoints} from "./symmetry_points";

export class ReciprocalLattice extends Lattice {

    constructor(config) {
        super(config);
    }

    get reciprocalVectors() {
        const vectors_ = this.vectors.vectorArrays;
        const divider = math.multiply(vectors_[0], math.cross(vectors_[1], vectors_[2]));
        return [
            math.multiply(math.cross(vectors_[1], vectors_[2]), 1 / divider),
            math.multiply(math.cross(vectors_[0], vectors_[2]), 1 / divider),
            math.multiply(math.cross(vectors_[1], vectors_[0]), 1 / divider)
        ];
    }

    getCartesianCoordinates(point) {
        return math.multiply(point, math.inv(this.reciprocalVectors))
    }

    get symmetryPoints() {
        return symmetryPoints(this);
    }

    get defaultKpointPath() {
        return paths[this.typeExtended] || paths[this.type];
    }

    extractKpointPath (dataPoints = []) {
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
