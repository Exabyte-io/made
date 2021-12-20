import math from "mathjs";

import { tolerance as TOLERANCE } from "./constants";

/**
 * @summary Zero threshold. Numbers below it are put to zero exactly.
 * Used to avoid math.js bug in treating zero as X.XXe-16.
 */
const EPSILON = 1e-8;
/**
 * @summary Returns scalar product of vectors
 * @param v1 {Number[]} Vector 1
 * @param v2 {Number[]} Vector 2
 * @return {Number}
 */

const product = (v1, v2) => {
    return math.multiply(v1, math.transpose(v2));
};

/**
 * @summary Returns length of a vector.
 * @param v {Number[]} Vector
 * @return {Number}
 */
const vlen = (v) => {
    return math.sqrt(product(v, v));
};

/**
 * @summary Returns angle between `a` and `b` vectors.
 * @param a {Number[]} Vector a
 * @param b {Number[]} Vector b
 * @param [unit] {String} `rad`, `deg`
 * @return {Number}
 */
const angle = (a, b, unit) => {
    const lenA = vlen(a);
    const lenB = vlen(b);
    return math.unit(math.acos(product(a, b) / (lenA * lenB)), "rad").toNumber(unit || "deg");
};

const angleUpTo90 = (...args) => {
    const angleUpTo180 = angle(...args);
    return angleUpTo180 < 90 ? angleUpTo180 : 180 - angleUpTo180;
};

/**
 * @summary Returns distance between 2 vectors.
 * @param v1 {Number[]} Vector
 * @param v2 {Number[]} Vector
 * @return {Number}
 */
const vDist = (v1, v2) => {
    if (v1.length !== v2.length) {
        console.error(
            "Attempting to calculate distance between vectors of different dimensionality",
        );
        return;
    }
    return vlen(v1.map((coordinate, index) => coordinate - v2[index]));
};

/**
 * @summary Returns checks whether 2 vector are equal within tolerance.
 * @param vec1 {Number[]} Vector
 * @param vec2 {Number[]} Vector
 * @param tolerance {Number} Tolerance
 * @return {Number}
 */
const vEqualWithTolerance = (vec1, vec2, tolerance = TOLERANCE.pointsDistance) =>
    vDist(vec1, vec2) <= tolerance;

/**
 * @summary Returns 0 if passed number is less than Made.math.EPSILON.
 * @param n {Number}
 * @return {Number}
 */
const roundToZero = (n) => {
    return Math.abs(n) < EPSILON ? 0 : n;
};

/**
 * @summary Returns number with specified precision.
 * @param x {Number}
 * @param n {Number}
 * @return {Number}
 */
const precise = (x, n = 7) => {
    return Number(x.toPrecision(n));
};

/**
 * @summary Returns mod of the passed value with the specified tolerance.
 * @param num {Number}
 * @param tolerance {Number}
 * @return {Number}
 */
const mod = (num, tolerance = 0.001) => {
    const m = num % 1;
    const x = num >= 0 ? m : 1 + m;

    if (math.smallerEq(Math.abs(x - 1), tolerance) || math.smallerEq(Math.abs(x), tolerance)) {
        return 0;
    }
    return x;
};

/**
 * @summary Returns cartesian of passed arrays.
 * @example combinations([1,2], [4,5], [6]) = [[1,4,6], [1,5,6], [2,4,6], [2,5,6]];
 */
const cartesianProduct = (...arg) => {
    const r = [];
    const max = arg.length - 1;

    const helper = (arr, i) => {
        for (let j = 0, l = arg[i].length; j < l; j++) {
            const a = arr.slice(0); // clone arr
            a.push(arg[i][j]);
            if (i === max) {
                r.push(a);
            } else {
                helper(a, i + 1);
            }
        }
    };

    helper([], 0);
    return r;
};

/**
 * @summary Returns all possible positive integer combinations where each value changes from 0 to a, b, c.
 * @param a {Number}
 * @param b {Number}
 * @param tolerance {Number}
 */
const almostEqual = (a, b, tolerance = TOLERANCE.pointsDistance) => {
    return Math.abs(a - b) < tolerance;
};

/**
 * @summary Returns true if number is 0 <= x < 1, inclusive, otherwise false.
 * Helper to deal with JS arithmetic artifacts.
 * @number number {Number}
 */
const isBetweenZeroInclusiveAndOne = (number, tolerance = TOLERANCE.length) => {
    return roundToZero(number) >= 0 && !almostEqual(number, 1, tolerance) && number < 1;
};

/**
 * @summary Returns all possible positive integer combinations where each value changes from 0 to a, b, c.
 * @example
 *   var comb = combinations(1, 2, 0);
 *   // [[0, 0, 0], [0, 1, 0], [0, 2, 0], [1, 0, 0], [1, 1, 0], [1, 2, 0]]
 * @param a
 * @param b
 * @param c
 */
const combinations = (a, b, c) => {
    const combs = [];
    for (let i = 0; i <= a; i++) {
        for (let j = 0; j <= b; j++) {
            for (let k = 0; k <= c; k++) {
                combs.push([i, j, k]);
            }
        }
    }
    return combs;
};

/*
 * @summary Same as `combinations` but accepting intervals (tuples) of integers: eg. [-3, 4]
 */
const combinationsFromIntervals = (arrA, arrB, arrC) => {
    const combs = [];
    for (let i = arrA[0]; i <= arrA[1]; i++) {
        for (let j = arrB[0]; j <= arrB[1]; j++) {
            for (let k = arrC[0]; k <= arrC[1]; k++) {
                combs.push([i, j, k]);
            }
        }
    }
    return combs;
};

export default {
    ...math,
    PI: Math.PI,
    trunc: Math.trunc,
    product,
    vlen,
    angle,
    angleUpTo90,
    vDist,
    vEqualWithTolerance,
    roundToZero,
    precise,
    mod,
    isBetweenZeroInclusiveAndOne,
    cartesianProduct,
    almostEqual,
    combinations,
    combinationsFromIntervals,
};
