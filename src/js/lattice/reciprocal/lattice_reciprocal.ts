import { ATOMIC_COORD_UNITS, units as UNITS } from "@mat3ra/code/dist/js/constants";
import { Vector3DSchema } from "@mat3ra/esse/dist/js/types";
import lodash from "lodash";

import math from "../../math";
import { Lattice } from "../lattice";
import { paths } from "./paths";
import { symmetryPoints } from "./symmetry_points";

export type KPointCoordinates = number[];
export type KPointPath = Array<{
    point: string;
    steps: number;
    coordinates: KPointCoordinates;
}>;

interface SymmetryPoint {
    point: string;
    coordinates: KPointCoordinates;
}

interface ConversionTable {
    [key: string]: {
        [key: string]: number;
    };
}

interface PathsType {
    [key: string]: Array<{
        point: string;
        steps: number;
    }>;
}

export class ReciprocalLattice extends Lattice {
    /**
     * Get reciprocal vectors for the current Lattice in cartesian (2pi / a) units
     * @return {Vector3DSchema[]}
     */
    get reciprocalVectors(): Vector3DSchema[] {
        const vectors_: Vector3DSchema[] = this.vectors.vectorArrays as Vector3DSchema[];
        const a: number = math.vlen(vectors_[0]);
        const divider: number =
            (math.multiply(
                vectors_[0],
                math.cross(vectors_[1], vectors_[2]),
            ) as unknown as number) / a;
        return [
            math.multiply(math.cross(vectors_[1], vectors_[2]), 1 / divider) as Vector3DSchema,
            math.multiply(math.cross(vectors_[2], vectors_[0]), 1 / divider) as Vector3DSchema,
            math.multiply(math.cross(vectors_[0], vectors_[1]), 1 / divider) as Vector3DSchema,
        ];
    }

    /**
     * Norms of reciprocal vectors.
     * @return {number[]}
     */
    get reciprocalVectorNorms(): number[] {
        return this.reciprocalVectors.map((vec) => math.norm(vec) as number);
    }

    /**
     * Ratio of reciprocal vector norms scaled by the inverse of the largest component.
     * @return {number[]}
     */
    get reciprocalVectorRatios(): number[] {
        const norms: number[] = this.reciprocalVectorNorms;
        const maxNorm: number = math.max(...norms) as number;
        return norms.map((n: number) => n / maxNorm);
    }

    /**
     * Get point (in crystal coordinates) in cartesian coordinates.
     * @param {KPointCoordinates} point - point in 3D space
     * @return {KPointCoordinates}
     */
    getCartesianCoordinates(point: KPointCoordinates): KPointCoordinates {
        return math.multiply(point, this.reciprocalVectors) as KPointCoordinates;
    }

    /**
     * Get the list of high-symmetry points for the current lattice.
     * @return {SymmetryPoint[]}
     */
    get symmetryPoints(): SymmetryPoint[] {
        return symmetryPoints(this);
    }

    /**
     * Get the default path in reciprocal space for the current lattice.
     * @return {Array<{point: string; steps: number}>}
     */
    get defaultKpointPath(): Array<{ point: string; steps: number }> {
        const pathsObj = paths as PathsType;
        return pathsObj[this.typeExtended as string] || pathsObj[this.type as string] || [];
    }

    /**
     * Find/mark the high symmetry points on a list with raw data and return the edited list.
     * @param {KPointCoordinates[]} dataPoints - list of point coordinates
     * @return {KPointPath}
     */
    extractKpointPath(dataPoints: KPointCoordinates[] = []): KPointPath {
        const kpointPath: KPointPath = [];
        const symmPoints: SymmetryPoint[] = this.symmetryPoints;

        dataPoints.forEach((point: KPointCoordinates, index: number) => {
            const symmPoint: SymmetryPoint | undefined = symmPoints.find((x) => {
                return math.vEqualWithTolerance(x.coordinates, point, 1e-4);
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
    calculateDimension(nPoints: number, index: number): number {
        const norms: number[] = this.reciprocalVectorNorms;
        const [j, k] = [0, 1, 2].filter((i) => i !== index); // get indices of other two dimensions
        const N: number = math.cbrt(
            (nPoints * norms[index] ** 2) / (norms[j] * norms[k]),
        ) as number;
        return math.max(1, math.ceil(N)) as number;
    }

    /**
     * Calculate grid dimensions from total number of k-points.
     * @param {number} nKpoints - Total number of k-points.
     * @return {number[]} - Grid dimensions
     */
    getDimensionsFromPointsCount(nKpoints: number): number[] {
        const indices: number[] = [0, 1, 2];
        return indices.map((i) => this.calculateDimension(nKpoints, i));
    }

    get conversionTable(): ConversionTable {
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
    getDimensionsFromSpacing(
        spacing: number,
        units: string = ATOMIC_COORD_UNITS.cartesian,
    ): number[] {
        const factor: number = this.conversionTable[units][ATOMIC_COORD_UNITS.cartesian] || 1;
        return this.reciprocalVectorNorms.map((norm: number) => {
            return math.max(1, math.ceil(lodash.round(norm / (spacing * factor), 4))) as number;
        });
    }

    /**
     * Calculate grid spacing as average of spacing along individual reciprocal axes.
     * @param {number[]} dimensions - Array of dimensions
     * @param {string} units - units of spacing parameter (default: 2pi / a)
     * @return {number} - average grid spacing
     */
    getSpacingFromDimensions(
        dimensions: number[],
        units: string = ATOMIC_COORD_UNITS.cartesian,
    ): number {
        const factor: number = this.conversionTable[ATOMIC_COORD_UNITS.cartesian][units] || 1;
        const norms: number[] = this.reciprocalVectorNorms;
        return (
            factor *
            (math.mean(
                dimensions.map((dim: number, i: number) => norms[i] / math.max(1, dim)),
            ) as number)
        );
    }
}
