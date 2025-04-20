import { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import { AtomicConstraintsSchema } from "@mat3ra/esse/dist/js/types";
import s from "underscore.string";

import { AtomicConstraints, AtomicConstraintValue } from "../constraints/constraints";
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

    constraints: AtomicConstraintsSchema;

    constructor(config: ConstrainedBasisConfig) {
        super(config);
        this.constraints = config.constraints;
        this._constraints = AtomicConstraints.fromObjects(config.constraints); // `constraints` is an Array with ids
    }

    static fromElementsCoordinatesAndConstraints(
        config: ElementsCoordinatesAndConstraintsConfig,
    ): ConstrainedBasis {
        const base = super.fromElementsAndCoordinates({
            elements: config.elements,
            coordinates: config.coordinates,
            units: config.units,
            cell: config.cell,
            labels: config.labels,
        }) as ConstrainedBasis;

        base.constraints = AtomicConstraints.fromValues(
            config.constraints,
        ).toJSON() as AtomicConstraintsSchema;
        base._constraints = AtomicConstraints.fromObjects(base.constraints);
        return base;
    }

    get AtomicConstraints() {
        return AtomicConstraints.fromObjects(this.constraints);
    }

    override toJSON(): AnyObject {
        return {
            ...super.toJSON(),
            constraints: this._constraints.toJSON(),
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
        return this.elementsCoordinatesConstraintsArray.map((entry) => {
            const element = entry[0] + entry[3]; // element with label, Fe1
            const coordinate = entry[1];
            const constraint = entry[2];
            return (
                s.sprintf("%-4s", element) +
                coordinate.prettyPrint() +
                " " +
                constraint.map((x) => (x ? 1 : 0)).join(" ")
            );
        });
    }
}
