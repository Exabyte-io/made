export const coefficients = {
    BOHR_TO_ANGSTROM: 0.52917721092,
    ANGSTROM_TO_BOHR: 1 / 0.52917721092
};

export const tolerance = {
    // in crystal coordinates
    length: 0.01,
    lengthAngstrom: 0.001,
    pointsDistance: 0.001
};

export const units = {
    bohr: 'bohr',
    angstrom: 'angstrom',
    degree: 'degree',
    radian: 'radian',
    alat: 'alat'
};

/**
 * @summary Coordinates units for a material's basis.
 */
export const ATOMIC_COORD_UNITS = {
    crystal: 'crystal',
    cartesian: 'cartesian'
};

// Only 3 digits will be considered for lattice and basis params on hashing
export const HASH_TOLERANCE = 3;


export default {
    coefficients,
    tolerance,
    units
}
