import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import type { AtomicConstraintsSchema, ConsistencyCheck, DerivedPropertiesSchema, FileSourceSchema, LatticeSchema, MaterialSchema } from "@mat3ra/esse/dist/js/types";
import type { BasisConfig } from "./basis/basis";
import { type ConstrainedBasisConfig, ConstrainedBasis } from "./basis/constrained_basis";
import { Constraint } from "./constraints/constraints";
import { Lattice } from "./lattice/lattice";
export declare const defaultMaterialConfig: MaterialSchema;
export interface MaterialSchemaJSON extends MaterialSchema, AnyObject {
}
export type MaterialMixinProps = ReturnType<typeof materialMixin>;
type MaterialMixinStaticProps = ReturnType<typeof materialMixinStaticProps>;
export type MaterialInMemoryEntity = InMemoryEntity & MaterialMixinProps;
export type MaterialMixinConstructor = Constructor<MaterialMixinProps> & MaterialMixinStaticProps;
export type OptionallyConstrainedBasisConfig = BasisConfig & Partial<Pick<ConstrainedBasisConfig, "constraints">>;
type Base = InMemoryEntity & NamedInMemoryEntity;
export declare function materialMixin<T extends Base = Base>(item: T): {
    toJSON(): MaterialSchema;
    name: string;
    src: FileSourceSchema | undefined;
    updateFormula(): void;
    /**
     * Gets Bolean value for whether or not a material is non-periodic vs periodic.
     * False = periodic, True = non-periodic
     */
    isNonPeriodic: boolean;
    /**
     * @summary Returns the specific derived property (as specified by name) for a material.
     */
    getDerivedPropertyByName(name: string): {
        name?: "volume" | undefined;
        units?: "angstrom^3" | undefined;
        value: number;
    } | {
        name?: "density" | undefined;
        units?: "g/cm^3" | undefined;
        value: number;
    } | {
        pointGroupSymbol?: string | undefined;
        spaceGroupSymbol?: string | undefined;
        tolerance?: {
            units?: "angstrom" | undefined;
            value: number;
        } | undefined;
        name?: "symmetry" | undefined;
    } | {
        name?: "elemental_ratio" | undefined;
        value: number;
        element?: string | undefined;
    } | {
        name?: "p-norm" | undefined;
        degree?: number | undefined;
        value: number;
    } | {
        name?: "inchi" | undefined;
        value: string;
    } | {
        name?: "inchi_key" | undefined;
        value: string;
    } | undefined;
    /**
     * @summary Returns the derived properties array for a material.
     */
    getDerivedProperties(): DerivedPropertiesSchema;
    /**
     * Gets material's formula
     */
    readonly formula: string;
    readonly unitCellFormula: string;
    unsetFileProps(): void;
    /**
     * @param textOrObject Basis text or JSON object.
     * @param format Format (xyz, etc.)
     * @param unitz crystal/cartesian
     */
    setBasis(textOrObject: string | BasisConfig, format?: string, unitz?: string): void;
    setBasisConstraints(constraints: Constraint[]): void;
    setBasisConstraintsFromArrayOfObjects(constraints: AtomicConstraintsSchema): void;
    readonly basis: OptionallyConstrainedBasisConfig;
    readonly Basis: ConstrainedBasis;
    /**
     * High-level access to unique elements from material instead of basis.
     */
    readonly uniqueElements: ("H" | "He" | "Li" | "Be" | "B" | "C" | "N" | "O" | "F" | "Ne" | "Na" | "Mg" | "Al" | "Si" | "P" | "S" | "Cl" | "Ar" | "K" | "Ca" | "Sc" | "Ti" | "V" | "Cr" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Zn" | "Ga" | "Ge" | "As" | "Se" | "Br" | "Kr" | "Rb" | "Sr" | "Y" | "Zr" | "Nb" | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" | "In" | "Sn" | "Sb" | "Te" | "I" | "Xe" | "Cs" | "Ba" | "La" | "Ce" | "Pr" | "Nd" | "Pm" | "Sm" | "Eu" | "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm" | "Yb" | "Lu" | "Hf" | "Ta" | "W" | "Re" | "Os" | "Ir" | "Pt" | "Au" | "Hg" | "Tl" | "Pb" | "Bi" | "Po" | "At" | "Rn" | "Fr" | "Ra" | "Ac" | "Th" | "Pa" | "U" | "Np" | "Pu" | "Am" | "Cm" | "Bk" | "Cf" | "Es" | "Fm" | "Md" | "No" | "Lr" | "Rf" | "Db" | "Sg" | "Bh" | "Hs" | "Mt" | "Ds" | "Rg" | "Cn" | "Nh" | "Fl" | "Mc" | "Lv" | "Ts" | "Og" | "X" | "Vac")[];
    lattice: LatticeSchema;
    readonly Lattice: Lattice;
    /**
     * Returns the inchi string from the derivedProperties for a non-periodic material, or throws an error if the
     *  inchi cannot be found.
     *  @returns {String}
     */
    getInchiStringForHash(): string;
    /**
     * Calculates hash from basis and lattice. Algorithm expects the following:
     * - asserts lattice units to be angstrom
     * - asserts basis units to be crystal
     * - asserts basis coordinates and lattice measurements are rounded to hash precision
     * - forms strings for lattice and basis
     * - creates MD5 hash from basisStr + latticeStr + salt
     * @param salt Salt for hashing, empty string by default.
     * @param isScaled Whether to scale the lattice parameter 'a' to 1.
     */
    calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
    hash: string;
    /**
     * Calculates hash from basis and lattice as above + scales lattice properties to make lattice.a = 1
     */
    readonly scaledHash: string;
    /**
     * Converts basis to crystal/fractional coordinates.
     */
    toCrystal(): void;
    /**
     * Converts current material's basis coordinates to cartesian.
     * No changes if coordinates already cartesian.
     */
    toCartesian(): void;
    /**
     * Returns material's basis in XYZ format.
     */
    getBasisAsXyz(fractional?: boolean): string;
    /**
     * Returns material in Quantum Espresso output format:
     * ```
     *    CELL_PARAMETERS (angstroms)
     *    -0.543131284  -0.000000000   0.543131284
     *    -0.000000000   0.543131284   0.543131284
     *    -0.543131284   0.543131284   0.000000000
     *
     *    ATOMIC_POSITIONS (crystal)
     *    Si       0.000000000   0.000000000  -0.000000000
     *    Si       0.250000000   0.250000000   0.250000000
     * ```
     */
    getAsQEFormat(): string;
    /**
     * Returns material in POSCAR format. Pass `true` to ignore original poscar source and re-serialize.
     */
    getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
    /**
     * Returns a copy of the material with conventional cell constructed instead of primitive.
     */
    getACopyWithConventionalCell(): T;
    /**
     * @summary a series of checks for the material and returns an array of results in ConsistencyChecks format.
     * @returns Array of checks results
     */
    getConsistencyChecks(): ConsistencyCheck[];
    /**
     * @summary a series of checks for the material's basis and returns an array of results in ConsistencyChecks format.
     * @returns Array of checks results
     */
    getBasisConsistencyChecks(): ConsistencyCheck[];
} & {
    _json: MaterialSchemaJSON;
};
export declare function materialMixinStaticProps<T extends Constructor<Base>>(item: T): {
    readonly defaultConfig: MaterialSchema;
    constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): FileSourceSchema;
};
export {};
