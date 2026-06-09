"use strict";
/**
 * Constants for molecule handling and precision.
 * Source of truth: /constants.json at repository root
 *
 * These values are duplicated here for TypeScript compilation.
 * When updating these, also update constants.json.
 */
Object.defineProperty(exports, "__esModule", { value: true });
exports.PRECISION_MAP = exports.molecularLatticePaddingFactor = exports.diatomicLatticePaddingFactor = exports.defaultNonPeriodicMinimumLatticeSize = void 0;
exports.defaultNonPeriodicMinimumLatticeSize = 3.0;
exports.diatomicLatticePaddingFactor = 3.0;
exports.molecularLatticePaddingFactor = 2.0;
exports.PRECISION_MAP = {
    coordinatesCrystal: 9,
    coordinatesCartesian: 6,
    latticeLengths: 6,
    latticeAngles: 4,
};
