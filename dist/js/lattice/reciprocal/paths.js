"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.paths = void 0;
/**
 * Default kpoing paths according to:
 *  [AFLOW](https://arxiv.org/abs/1004.2974) methodology.
 *  Paths are split in parts for clarity.
 */
const points = {
    CUB: [
        ["Г", "X", "M", "Г", "R", "X"],
        ["M", "R"],
    ],
    BCC: [
        ["Г", "H", "N", "Г", "P", "H"],
        ["P", "N"],
    ],
    FCC: [
        ["Г", "X", "W", "K", "Г", "L", "U", "W", "L"],
        ["U", "X"],
    ],
    TET: [
        ["Г", "X", "M", "Г", "Z", "R", "A", "Z"],
        ["X", "R"],
        ["M", "A"],
    ],
    "BCT-1": [
        ["Г", "X", "M", "Г", "Z", "P", "N", "Z1", "M"],
        ["X", "P"],
    ],
    "BCT-2": [
        ["Г", "X", "Y", "∑", "Г", "Z", "∑1", "N", "P", "Y1", "Z"],
        ["X", "P"],
    ],
    ORC: [
        ["Г", "X", "S", "Y", "Г", "Z", "U", "R", "T", "Z"],
        ["Y", "T"],
        ["U", "X"],
        ["S", "R"],
    ],
    "ORCF-1": [
        ["Г", "Y", "T", "Z", "Г", "X", "A1", "Y"],
        ["T", "X1"],
        ["X", "A", "Z"],
        ["L", "Г"],
    ],
    "ORCF-2": [
        ["Г", "Y", "C", "D", "X", "Г", "Z", "D1", "H", "C"],
        ["C1", "Z"],
        ["X", "H1"],
        ["H", "Y"],
        ["L", "Г"],
    ],
    "ORCF-3": [
        ["Г", "Y", "T", "Z", "Г", "X", "A1", "Y"],
        ["X", "A", "Z"],
        ["L", "Г"],
    ],
    ORCI: [
        ["Г", "X", "L", "T", "W", "R", "X1", "Z", "Г", "Y", "S", "W"],
        ["L1", "Y"],
        ["Y1", "Z"],
    ],
    ORCC: [
        ["Г", "X", "S", "R", "A", "Z", "Г", "Y", "X1", "A1", "T", "Y"],
        ["Z", "T"],
    ],
    HEX: [
        ["Г", "M", "K", "Г", "A", "L", "H", "A"],
        ["L", "M"],
        ["K", "H"],
    ],
    "RHL-1": [
        ["Г", "L", "B1"],
        ["B", "Z", "Г", "X"],
        ["Q", "F", "P1", "Z"],
        ["L", "P"],
    ],
    "RHL-2": [["Г", "P", "Z", "Q", "Г", "F", "P1", "Q1", "L", "Z"]],
    MCL: [
        ["Г", "Y", "H", "C", "E", "M1", "A", "X", "H1"],
        ["M", "D", "Z"],
        ["Y", "D"],
    ],
    "MCLC-1": [
        ["Г", "Y", "Г", "L", "I"],
        ["I1", "Z", "F1"],
        ["Y", "X1"],
        ["X", "Г", "N"],
        ["M", "Г"],
    ],
    "MCLC-2": [
        ["Г", "Y", "F", "L", "I"],
        ["I1", "Z", "F1"],
        ["N", "Г", "M"],
    ],
    "MCLC-3": [
        ["Г", "Y", "F", "H", "Z", "I", "F1"],
        ["H1", "Y1", "X", "Г", "N"],
        ["M", "Г"],
    ],
    "MCLC-4": [
        ["Г", "Y", "F", "H", "Z", "I"],
        ["H1", "Y1", "X", "Г", "N"],
        ["M", "Г"],
    ],
    "MCLC-5": [
        ["Г", "Y", "F", "L", "I"],
        ["I1", "Z", "H", "F1"],
        ["H1", "Y1", "X", "Г", "N"],
        ["M", "Г"],
    ],
    TRI: [
        ["X", "Г", "Y"],
        ["L", "Г", "Z"],
        ["N", "Г", "M"],
        ["R", "Г"],
    ],
};
// Export a processed version of the paths
exports.paths = (() => {
    const result = {};
    Object.entries(points).forEach(([key, pathSegments]) => {
        // Flatten arrays of path segments into a single array
        const flattenedPath = pathSegments.reduce((acc, segment) => [...acc, ...segment], []);
        // Convert to KPointStep array
        // TODO: calculate number of steps based on distance in k-space
        result[key] = flattenedPath.map((point) => ({
            point,
            steps: 10,
        }));
    });
    return result;
})();
