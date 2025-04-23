import { AtomicConstraintsSchema } from "@mat3ra/esse/dist/js/types";

import { AtomicConstraints, AtomicConstraintValue, Constraint } from "../constraints/constraints";
import { Basis, BasisConfig, ElementsAndCoordinatesConfig } from "./basis";
import { Coordinate } from "./coordinates";

export interface ConstrainedBasisConfig extends BasisConfig {
    constraints: AtomicConstraintsSchema;
}

export interface ElementsCoordinatesAndConstraintsConfig extends ElementsAndCoordinatesConfig {
    constraints: AtomicConstraintValue[];
}

/**
 * @summary Extension of the Basis class able to deal with atomic constraints.
 * @extends Basis
 */
export class ConstrainedBasis extends Basis {
    _constraints: AtomicConstraints;

    constructor(config: ConstrainedBasisConfig) {
        super(config);
        const { constraints } = config;
        this._constraints = AtomicConstraints.fromObjects(constraints || []); // `constraints` is an Array with ids
    }

    static fromElementsCoordinatesAndConstraints(
        config: ElementsCoordinatesAndConstraintsConfig,
    ): ConstrainedBasis {
        const basisConfig = this._convertValuesToConfig(config);
        const constraints = AtomicConstraints.fromValues(config.constraints);
        return new this({
            ...basisConfig,
            constraints: constraints.toJSON() as AtomicConstraintsSchema,
        });
    }

    get constraints() {
        return this._constraints.toJSON() as AtomicConstraintsSchema;
    }

    set constraints(constraints: AtomicConstraintsSchema) {
        this._constraints = AtomicConstraints.fromObjects(constraints);
    }

    get AtomicConstraints() {
        return AtomicConstraints.fromObjects(this.constraints);
    }

    override toJSON(): ConstrainedBasisConfig {
        return {
            ...super.toJSON(),
            constraints: this.constraints,
        };
    }

    getConstraintByIndex(idx: number): AtomicConstraintValue {
        return this._constraints.getElementValueByIndex(idx) || [false, false, false];
    }

    /**
     * Helper function returning a nested array with [element, coordinates, constraints] as elements
     */
    get elementsCoordinatesConstraintsArray(): [
        string,
        Coordinate,
        AtomicConstraintValue,
        string,
    ][] {
        return this._elements.values.map((element: any, idx: number) => {
            const coordinate = this.getCoordinateByIndex(idx);
            const constraint = this.getConstraintByIndex(idx);
            const atomicLabel = this.atomicLabelsArray[idx];
            return [element, coordinate, constraint, atomicLabel];
        });
    }

    /**
     * Returns an array with atomic positions (with constraints) per atom stored as strings.
     * E.g., ``` ['Si  0 0 0  0 1 0', 'Li  0.5 0.5 0.5  1 0 1']```
     */
    get atomicPositionsWithConstraints(): string[] {
        return this.elementsCoordinatesConstraintsArray.map(
            ([element, coordinate, constraint, label]) => {
                const fullElement = element + label; // e.g., Fe1
                return new Constraint({ id: 0, value: constraint }).prettyPrint(
                    fullElement,
                    coordinate,
                    constraint,
                );
            },
        );
    }
}
