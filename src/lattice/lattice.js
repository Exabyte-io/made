import lodash from "lodash";
import math from "../math";
import constants, {HASH_TOLERANCE} from "../constants";

import {LATTICE_TYPE, LATTICE_TYPE_EXTENDED, LATTICE_TYPE_CONFIGS} from "./types";
import {LatticeBravais} from "./lattice_bravais";
import {LatticeVectors} from "./lattice_vectors";
import {primitiveCell} from "../cell/primitive_cell";
import {Cell} from "../cell/cell";

/*
 * @summary Container class for crystal lattice and associated methods.
 *          Follows Bravais convention for lattice types and contains lattice vectors within.
 *          Units for lattice vector coordinates are "angstroms", and "degrees" for the corresponding angles.
 */
export class Lattice extends LatticeBravais {

    constructor(config = {}) {
        super(config);
        this.vectors = LatticeVectors.fromBravais(config);
    }

    static fromVectors(config) {
        return new Lattice(LatticeBravais.fromVectors(config).toJSON());
    }

    toJSON(skipRounding = false) {
        const round = skipRounding ? () => {} : Lattice._roundValue;  // round values by default
        return {
            a: round(this.a),
            b: round(this.b),
            c: round(this.c),
            alpha: round(this.alpha),
            beta: round(this.beta),
            gamma: round(this.gamma),
            units: {
                length: this.units.length,
                angle: this.units.angle,
            },
            type: this.type,
            vectors: this.vectors.toJSON(),
        }
    }

    clone(extraContext) {return new this.constructor(Object.assign({}, this.toJSON(), extraContext))}

    get vectorArrays() {
        return this.vectors.vectorArrays;
    }

    get Cell() {
        return new Cell(this.vectorArrays);
    }

    get typeLabel() {
        return lodash.get(LATTICE_TYPE_CONFIGS.find(c => c.code === this.type), 'label', 'Unknown');
    }

    get typeExtended() {
        const {a, b, c, alpha, beta, gamma, type} = this;
        let cosAlpha = math.cos(alpha / 180 * math.PI);

        switch (type) {
            case LATTICE_TYPE.BCT:
                return c < a ? LATTICE_TYPE_EXTENDED.BCT_1 : LATTICE_TYPE_EXTENDED.BCT_2;
            case LATTICE_TYPE.ORCF:
                if (1 / (a * a) >= 1 / (b * b) + 1 / (c * c)) {
                    return LATTICE_TYPE_EXTENDED.ORCF_1;
                }
                return LATTICE_TYPE_EXTENDED.ORCF_2;
            case LATTICE_TYPE.RHL:
                return cosAlpha > 0 ? LATTICE_TYPE_EXTENDED.RHL_1 : LATTICE_TYPE_EXTENDED.RHL_2;
            case LATTICE_TYPE.MCLC:
                if (gamma >= 90) {
                    // MCLC-1,2
                    return LATTICE_TYPE_EXTENDED.MCLC_1;
                } else if (b / c * cosAlpha + b * b / (a * a) * (1 - cosAlpha * cosAlpha) <= 1) {
                    // MCLC-3,4
                    return LATTICE_TYPE_EXTENDED.MCLC_3;
                }
                return LATTICE_TYPE_EXTENDED.MCLC_5;
            case LATTICE_TYPE.TRI:
                if (alpha > 90 && beta > 90 && gamma >= 90) {
                    // TRI-1a,2a
                    return LATTICE_TYPE_EXTENDED.TRI_1a;
                }
                return LATTICE_TYPE_EXTENDED.TRI_1b;
            default:
                return type;
        }
    }

    get volume() {
        return math.abs(math.det(this.vectors));
    }

    /*
     * @summary: Returns a "default" primitive lattice for each lattice type with lattice parameters scaled by "a",
     * @param latticeConfig {Object} LatticeBravais config (see constructor)
     */
    static getDefaultPrimitiveLatticeConfigByType(latticeConfig) {
        const f_ = Lattice._roundValue;
        // construct new primitive cell using lattice parameters and skip rounding the vectors
        const pCell = primitiveCell(latticeConfig, true);
        // create new lattice from primitive cell
        const newLattice = Object.assign({}, Lattice.fromVectorArrays(pCell, latticeConfig.type));

        // preserve the new type and scale back the lattice parameters
        let k = latticeConfig.a / newLattice.a;

        return Object.assign(newLattice, {
            a: f_(newLattice.a * k),
            b: f_(newLattice.b * k),
            c: f_(newLattice.c * k),
            alpha: f_(newLattice.alpha),
            beta: f_(newLattice.beta),
            gamma: f_(newLattice.gamma),
        });

    }

    get unitCell() {
        const vectors = _.flatten(this.vectorArrays);
        vectors.push(this.units.length);
        return new UnitCell(...vectors);

    }

    // returns a string for hash calculation
    getHashString(isScaled = false) {
        // lattice vectors must be measured in angstroms
        const latticeInAngstroms = this;
        let scaleK = (isScaled) ? latticeInAngstroms.a : 1;
        const scaledLattice = Object.assign({}, latticeInAngstroms, {
            a: latticeInAngstroms.a / scaleK,
            b: latticeInAngstroms.b / scaleK,
            c: latticeInAngstroms.c / scaleK
        });
        // form lattice string
        return [
            scaledLattice.a,
            scaledLattice.b,
            scaledLattice.c,
            scaledLattice.alpha,
            scaledLattice.beta,
            scaledLattice.gamma
        ].map(x => math.round(x, HASH_TOLERANCE)).join(';') + ";";
    }

}

// TODO: refactor and remove the need for UnitCell
class UnitCell {
    constructor(ax, ay, az, bx, by, bz, cx, cy, cz, units) {
        this.ax = ax;
        this.ay = ay;
        this.az = az;
        this.bx = bx;
        this.by = by;
        this.bz = bz;
        this.cx = cx;
        this.cy = cy;
        this.cz = cz;
        this.units = units;
    }

    vectorA() {
        return [this.ax, this.ay, this.az];
    }

    vectorB() {
        return [this.bx, this.by, this.bz];
    }

    vectorC() {
        return [this.cx, this.cy, this.cz];
    }

    axes() {
        return [this.vectorA(), this.vectorB(), this.vectorC()]
    }
}
