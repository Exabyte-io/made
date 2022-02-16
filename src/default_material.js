import { ATOMIC_COORD_UNITS, units } from "./constants";
import { LATTICE_TYPE } from "./lattice/types";

export const defaultMaterialConfig = {
    name: "Silicon FCC",
    basis: {
        elements: [
            {
                id: 1,
                value: "Si",
            },
            {
                id: 2,
                value: "Si",
            },
        ],
        coordinates: [
            {
                id: 1,
                value: [0.0, 0.0, 0.0],
            },
            {
                id: 2,
                value: [0.25, 0.25, 0.25],
            },
        ],
        units: ATOMIC_COORD_UNITS.crystal,
    },
    lattice: {
        // Primitive cell for Diamond FCC Silicon at ambient conditions
        type: LATTICE_TYPE.FCC,
        a: 3.867,
        b: 3.867,
        c: 3.867,
        alpha: 60,
        beta: 60,
        gamma: 60,
        units: {
            length: units.angstrom,
            angle: units.degree,
        },
    },
};

/*
 * Function returns the defaultMaterial object.
 * @returns {Object}
 */
export function getDefaultMaterialConfig() {
    return defaultMaterialConfig;
}

export default {
    defaultMaterialConfig,
    getDefaultMaterialConfig,
};
