import {
    Coordinate3DSchema,
    MaterialSchema,
    Matrix3X3Schema,
    Vector3DSchema,
} from "@mat3ra/esse/dist/js/types";

import { Cell } from "../cell/cell";
import { Lattice } from "../lattice/lattice";
import { Material } from "../material";
import math from "../math";
import SupercellTools from "./supercell";

const MULT = math.multiply;
const ADD = math.add;
const DOT = math.product;
const getMatrixInLeftHandedRepresentation = (matrix: Matrix3X3Schema) => {
    return math.det(matrix) < 0 ? MULT(matrix, -1) : matrix;
};

/**
 * Helper function for extended GCD.
 * Inspired by https://gitlab.com/ase/ase/-/blob/20401e2596821a85f17ced1cb18f32c508e35fb0/ase/build/general_surface.py
 * @param {Number} a
 * @param {Number} b
 * @return {Array}
 */
function extGCD(a: number, b: number): [number, number] {
    if (b === 0) return [1, 0];
    if (math.mod(a, b) === 0) return [0, 1];

    const [x, y] = extGCD(b, math.mod(a, b));
    return [y, x - y * math.floor(a / b)];
}

/**
 * Generates a slab scaling matrix for the specified cell based on miller indices.
 * Inspired by from https://gitlab.com/ase/ase/blob/master/ase/build/general_surface.py.
 * @param cell {Cell}
 * @param millerIndices {Number[]}
 * @param tol {Number} Zero-value tolerance
 * @return {Number[][]}
 */
function getMillerScalingMatrix(
    cell: Cell,
    millerIndices: Coordinate3DSchema,
    tol = 1e-8,
): Matrix3X3Schema {
    if (!millerIndices.reduce((a, b) => math.abs(a) + math.abs(b)))
        throw new Error("Miller indices are zeros.");

    let scalingMatrix;

    const [h, k, l] = millerIndices;
    const [h0, k0, l0] = millerIndices.map((i) => i === 0);

    if ((h0 && k0) || (h0 && l0) || (k0 && l0)) {
        // if any two indices are zero

        if (!h0) {
            scalingMatrix = [
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 0],
            ];
        }
        if (!k0) {
            scalingMatrix = [
                [0, 0, 1],
                [1, 0, 0],
                [0, 1, 0],
            ];
        }
        if (!l0) {
            scalingMatrix = [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ];
        }
    } else {
        const [a1, a2, a3] = cell.vectorArrays;
        // z1 = (k * a1 - h * a2)
        // z2 = (l * a1 - h * a3)
        // z3 = l * a2 - k * a3)
        const z1 = ADD(MULT(k, a1), -MULT(h, a2));
        const z2 = ADD(MULT(l, a1), -MULT(h, a3));
        const z3 = ADD(MULT(l, a2), -MULT(k, a3));

        let [p, q] = extGCD(k, l);

        // @ts-ignore
        const k1 = DOT(ADD(MULT(p, z1), MULT(q, z2)), z3);
        // @ts-ignore
        const k2 = DOT(ADD(MULT(l, z1), -MULT(k, z2)), z3);

        if (math.abs(k2) > tol) {
            // For mathjs version 3.20: round(-0.5) = -0
            // For mathjs version 5.10: round(-0.5) = -1
            // Here we specify rounding method to Bankers
            // For Python 3.11: round(-0.5) = 0
            const value = k1 / k2;
            const roundedValue = math.roundCustom(value, 0, math.RoundingMethod.Bankers);
            const i = -roundedValue;
            [p, q] = [p + i * l, q - i * k];
        }

        const [a, b] = extGCD(p * k + q * l, h);
        const c1 = [p * k + q * l, -p * h, -q * h];
        const c2 = [0, l, -k].map((c) => math.trunc(c / math.gcd(l, k))); // floor division
        const c3 = [b, a * p, a * q];

        scalingMatrix = [c1, c2, c3];
    }

    // @ts-ignore
    return getMatrixInLeftHandedRepresentation(scalingMatrix);
}

/**
 * Scale the surface cell in out-of-plane direction according to thickness
 *  and construct the lateral supercell from vx/vy parameters.
 * @param bulkCell {Cell}
 * @param surfaceCell {Cell}
 * @param outOfPlaneAxisIndex {Number} Index of the cell vector most collinear with the out-of-surface-plane vector
 * @param thickness {Number} Surface (Slab) thickness in layers (Positive Integer).
 * @param vx {Number} Size of lateral supercell along the direction of the first (x) cell vector (Positive Integer).
 * @param vy {Number} Size of lateral supercell along the direction of the second (y) cell vector (Positive Integer).
 * @return {Number[][]}
 */
function getDimensionsScalingMatrix(
    outOfPlaneAxisIndex: number,
    thickness: number,
    vx: number,
    vy: number,
): Matrix3X3Schema {
    const transformationMatrix = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ] as Matrix3X3Schema;
    const vxIndex = outOfPlaneAxisIndex === 2 ? 0 : outOfPlaneAxisIndex + 1;
    const vyIndex = vxIndex === 2 ? 0 : vxIndex + 1;

    transformationMatrix[outOfPlaneAxisIndex] = MULT(
        thickness,
        transformationMatrix[outOfPlaneAxisIndex],
    ) as Vector3DSchema;
    transformationMatrix[vxIndex] = MULT(vx, transformationMatrix[vxIndex]) as Vector3DSchema;
    transformationMatrix[vyIndex] = MULT(vy, transformationMatrix[vyIndex]) as Vector3DSchema;

    return transformationMatrix as Matrix3X3Schema;
}

export type SlabConfigSchema = MaterialSchema & {
    outOfPlaneAxisIndex: number;
};

/*
 * Generate a config object for a surface based on a set of parameters.
 * @param material {Material} The parent material (bulk) this surface is constructed from
 * @param millerIndices {Number[]} Miller Indices that define the surface cleavage.
 * @param numberOfLayers {Number} Surface (Slab) thickness in layers (Positive Integer).
 * @param vx {Number} Size of lateral supercell along the direction of the first (x) cell vector (Positive Integer).
 * @param vy {Number} Size of lateral supercell along the direction of the second (y) cell vector (Positive Integer).
 * @return {Object}
 */
function generateConfig(
    material: Material,
    millerIndices: Coordinate3DSchema,
    numberOfLayers = 1,
    vx = 1,
    vy = 1,
): SlabConfigSchema {
    if (numberOfLayers < 1)
        throw new Error("Made.tools.surface.generateConfig: number of layers < 1.");

    const cell = material.Lattice.vectors;
    const millerScalingMatrix = getMillerScalingMatrix(cell, millerIndices);
    const millerSupercell = cell.cloneAndScaleByMatrix(millerScalingMatrix);
    const millerPlanePseudoNormal = cell.convertPointToCartesian(millerIndices);
    const outOfPlaneAxisIndex =
        millerSupercell.getMostCollinearVectorIndex(millerPlanePseudoNormal);
    const dimensionsScalingMatrix = getDimensionsScalingMatrix(
        outOfPlaneAxisIndex,
        numberOfLayers,
        vx,
        vy,
    );
    const supercellMatrix = MULT(dimensionsScalingMatrix, millerScalingMatrix);
    const supercell = millerSupercell.cloneAndScaleByMatrix(dimensionsScalingMatrix);
    const tempBasis = material.Basis.clone();
    const newBasis = SupercellTools.generateNewBasisWithinSupercell(
        tempBasis,
        cell,
        supercell,
        supercellMatrix as Matrix3X3Schema,
    );
    const newLattice = Lattice.fromVectors({
        a: supercell.vectorArrays[0],
        b: supercell.vectorArrays[1],
        c: supercell.vectorArrays[2],
    });

    return {
        name: `${material.name} - slab ${JSON.stringify(millerIndices)}`,
        basis: newBasis.toJSON(),
        lattice: newLattice.toJSON(),
        // extra parameter for use in creating vacuum etc.
        outOfPlaneAxisIndex,
    };
}

export default {
    generateConfig,
};
