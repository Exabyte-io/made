/**
 * Constants for molecule handling and precision.
 * Source of truth: /constants.json at repository root
 *
 * These values are duplicated here for TypeScript compilation.
 * When updating these, also update constants.json.
 */
export declare const defaultNonPeriodicMinimumLatticeSize = 3;
export declare const diatomicLatticePaddingFactor = 3;
export declare const molecularLatticePaddingFactor = 2;
export declare const PRECISION_MAP: {
    coordinatesCrystal: number;
    coordinatesCartesian: number;
    latticeLengths: number;
    latticeAngles: number;
};
