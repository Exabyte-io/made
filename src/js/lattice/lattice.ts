import { HASH_TOLERANCE } from "@mat3ra/code/dist/js/constants";
import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import {
    Coordinate3DSchema,
    LatticeSchema,
    LatticeTypeEnum,
    LatticeTypeExtendedEnum,
    LatticeVectorsSchema,
    Matrix3X3Schema,
} from "@mat3ra/esse/dist/js/types";
import * as lodash from "lodash";

import { Cell } from "../cell/cell";
import { getPrimitiveLatticeVectorsFromConfig } from "../cell/primitive_cell";
import math from "../math";
import { LATTICE_TYPE_CONFIGS } from "./lattice_types";
import { UnitCell, UnitCellProps } from "./unit_cell";

/**
 * Scaling factor used to calculate the new lattice size for non-periodic systems.
 * The scaling factor ensures that a non-periodic structure will have have a lattice greater than the structures size.
 */
export const nonPeriodicLatticeScalingFactor = 2.0;

export class LatticeVectors extends Cell implements LatticeVectorsSchema {}

export type { LatticeVectorsSchema };

export class Lattice extends InMemoryEntity implements LatticeSchema {
    static defaultConfig: LatticeSchema = {
        a: 1,
        b: 1,
        c: 1,
        alpha: 90,
        beta: 90,
        gamma: 90,
        units: { length: "angstrom", angle: "degree" },
        type: "CUB",
    };

    a: LatticeSchema["a"];

    b: LatticeSchema["b"];

    c: LatticeSchema["c"];

    alpha: LatticeSchema["alpha"];

    beta: LatticeSchema["beta"];

    gamma: LatticeSchema["gamma"];

    type: LatticeTypeEnum = "TRI";

    units: LatticeSchema["units"];

    constructor(config: LatticeSchema = Lattice.defaultConfig) {
        super(config);
        const { a, b, c, alpha, beta, gamma, units, type } = config;
        this.a = a;
        this.b = b;
        this.c = c;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.units = units;
        this.type = type || "TRI";
    }

    static fromConfig(config: object): Lattice {
        const latticeConfig = config as LatticeSchema;
        return new Lattice(latticeConfig);
    }

    static fromConfigPartial(config: LatticeSchema): Lattice {
        const primitiveLatticeConfig = Lattice.getDefaultPrimitiveLatticeConfigByType(config);
        return new Lattice(primitiveLatticeConfig);
    }

    calculateVectors(): Matrix3X3Schema {
        const { a } = this;
        const { b } = this;
        const { c } = this;

        const alphaRad = math.unit(this.alpha, "deg").toNumber("rad");
        const betaRad = math.unit(this.beta, "deg").toNumber("rad");
        const gammaRad = math.unit(this.gamma, "deg").toNumber("rad");

        const cosAlpha = math.cos(alphaRad);
        const cosBeta = math.cos(betaRad);
        const cosGamma = math.cos(gammaRad);
        const sinAlpha = math.sin(alphaRad);
        const sinBeta = math.sin(betaRad);

        const gammaStar = math.acos((cosAlpha * cosBeta - cosGamma) / (sinAlpha * sinBeta));
        const cosGammaStar = math.cos(gammaStar);
        const sinGammaStar = math.sin(gammaStar);

        const vectorA: Coordinate3DSchema = [a * sinBeta, 0.0, a * cosBeta];
        const vectorB: Coordinate3DSchema = [
            -b * sinAlpha * cosGammaStar,
            b * sinAlpha * sinGammaStar,
            b * cosAlpha,
        ];
        const vectorC: Coordinate3DSchema = [0.0, 0.0, c];

        return [vectorA, vectorB, vectorC];
    }

    static fromVectors(config: LatticeVectorsSchema): Lattice {
        return this.fromVectorsArray([config.a, config.b, config.c]);
    }

    static fromVectorsArray(
        vectors: Matrix3X3Schema,
        units: LatticeSchema["units"] = Lattice.defaultConfig.units,
        type: LatticeTypeEnum = "TRI",
    ): Lattice {
        const [aVec, bVec, cVec] = vectors;
        const a = math.vlen(aVec);
        const b = math.vlen(bVec);
        const c = math.vlen(cVec);
        const alpha = math.angle(bVec, cVec, "deg");
        const beta = math.angle(aVec, cVec, "deg");
        const gamma = math.angle(aVec, bVec, "deg");

        return new Lattice({
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            units,
            type,
        });
    }

    // TODO: remove
    get unitCell() {
        const units = this.units?.length || "angstrom";
        const vectors = [...lodash.flatten(this.vectorArrays), units] as UnitCellProps;
        return new UnitCell(vectors);
    }

    get vectors(): LatticeVectors {
        return LatticeVectors.fromVectorsArray(this.calculateVectors());
    }

    get vectorArrays(): Matrix3X3Schema {
        return this.vectors.vectorArrays;
    }

    get vectorArraysRounded(): Matrix3X3Schema {
        return this.vectors.vectorArraysRounded;
    }

    get cellVolume(): number {
        return this.vectors.volume;
    }

    get cellVolumeRounded(): number {
        return this.vectors.volumeRounded;
    }

    /**
     * Get a short label for the type of the lattice, eg. "MCLC".
     */
    get typeLabel(): string {
        return lodash.get(
            LATTICE_TYPE_CONFIGS.find((c) => c.code === this.type),
            "label",
            "Unknown",
        );
    }

    /**
     * Get a short label for the extended type of the lattice, eg. "MCLC-5".
     */
    get typeExtended(): LatticeTypeExtendedEnum {
        const { a, b, c, alpha, beta, gamma, type } = this;
        const cosAlpha = math.cos((alpha / 180) * math.PI);

        switch (type) {
            case "BCT":
                return c < a ? "BCT-1" : "BCT-2";
            case "ORCF":
                if (1 / (a * a) >= 1 / (b * b) + 1 / (c * c)) {
                    return "ORCF-1";
                }
                return "ORCF-2";
            case "RHL":
                return cosAlpha > 0 ? "RHL-1" : "RHL-2";
            case "MCLC":
                if (gamma >= 90) {
                    // MCLC-1,2
                    return "MCLC-1";
                }
                if ((b / c) * cosAlpha + ((b * b) / (a * a)) * (1 - cosAlpha * cosAlpha) <= 1) {
                    // MCLC-3,4
                    return "MCLC-3";
                }
                return "MCLC-5";
            case "TRI":
                if (alpha > 90 && beta > 90 && gamma >= 90) {
                    // TRI-1a,2a
                    return "TRI_1a";
                }
                return "TRI_1b";
            default:
                return type;
        }
    }

    /**
     * Calculate the volume of the lattice cell.
     */
    get volume(): number {
        return math.abs(math.det(this.vectorArrays));
    }

    /*
     * Returns a "default" primitive lattice by type, with lattice parameters scaled by the length of "a",
     * @param latticeConfig {Object} LatticeBravais config (see constructor)
     */
    static getDefaultPrimitiveLatticeConfigByType(latticeConfig: LatticeSchema) {
        const f_ = math.roundArrayOrNumber;
        // construct new primitive cell using lattice parameters and skip rounding the vectors
        const vectors = getPrimitiveLatticeVectorsFromConfig(latticeConfig);
        // create new lattice from primitive cell
        const newLattice = Lattice.fromVectorsArray(vectors, undefined, latticeConfig.type);
        // preserve the new type and scale back the lattice parameters
        const k = latticeConfig.a / newLattice.a;

        return Object.assign(newLattice, {
            a: f_(newLattice.a * k),
            b: f_(newLattice.b * k),
            c: f_(newLattice.c * k),
            alpha: f_(newLattice.alpha),
            beta: f_(newLattice.beta),
            gamma: f_(newLattice.gamma),
        });
    }

    /**
     * Returns a string further used for the calculation of an unique hash.
     * @param isScaled - Whether to scale the vectors by the length of the first vector initially.
     */
    getHashString(isScaled = false): string {
        // lattice vectors must be measured in angstroms
        const scaleK = isScaled ? this.a : 1;
        const scaledLattice = {
            ...this,
            a: this.a / scaleK,
            b: this.b / scaleK,
            c: this.c / scaleK,
        };
        // form lattice string
        return `${[
            scaledLattice.a,
            scaledLattice.b,
            scaledLattice.c,
            scaledLattice.alpha,
            scaledLattice.beta,
            scaledLattice.gamma,
        ]
            .map((x) => math.round(x, HASH_TOLERANCE))
            .join(";")};`;
    }

    /**
     * Get the list of editable keys (eg. 'a', 'alpha') for the current lattice.
     * @return {Object}
     * @example {a: true, b: false, c: false, alpha: true, beta: false, gamma: false}
     */
    get editables() {
        const object = {};
        const editablesList = LATTICE_TYPE_CONFIGS.find(
            (entry: any) => entry.code === this.type,
        )?.editables;
        // ["a", "gamma"] => {a: true, gamma: true}
        if (editablesList) {
            editablesList.forEach((element: any) => {
                Object.assign(object, {
                    [element]: true,
                });
            });
        }

        return object;
    }

    // @ts-ignore
    toJSON(exclude?: string[]): LatticeSchema {
        return {
            ...super.toJSON(exclude),
            vectors: this.vectors.toJSON(),
        } as LatticeSchema;
    }
}
