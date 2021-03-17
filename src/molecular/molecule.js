import _ from "underscore";
import lodash from "lodash";

import math from "../math";
import constants, {HASH_TOLERANCE} from "../constants";

/**
 * Container class for molecular structures and associated methods.
 * Units for molecular structures in zmatrix are angstroms for bonds, degrees for angles & dihedrals
 */
export class Molecule {
    /**
     * Returns a string further used for the calculation of an unique hash.
     * @param {Boolean} isScaled - Whether to scale the vectors by the length of the first vector initially.
     * @return {string}
     */
    getZmatrixHashString(isScaled = false, zmatrix = []) {
        // Bonds need to be measured in angstroms, angles & dihedrals in degrees
        var bonds = zmatrix[0]
        var angles = zmatrix[1]
        var dihedrals = zmatrix[2]
        var sum_bonds = bonds.reduce(function(a, b ){
            return a + b;
        }, 0);
        var sum_angles = angles.reduce(function(a, b ){
            return a + b;
        }, 0);
        var sum_dihedrals = dihedrals.reduce(function(a, b ){
            return a + b;
        }, 0);

        const sumZmatrix = Object.assign({}, sum_bonds, sum_angles, sum_dihedrals, {
            a: sum_bonds,
            b: sum_angles,
            d: sum_dihedrals
        });
        return [
            sumZmatrix.a,
            sumZmatrix.b,
            sumZmatrix.d
        ].map(x => math.round(x, HASH_TOLERANCE)).join(';') + ";";
    }
}




