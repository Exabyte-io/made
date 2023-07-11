import { ConstrainedBasis } from "../basis/constrained_basis";
import { Lattice } from "../lattice/lattice";
import { BaseParser } from "./init";

export class MaterialParser extends BaseParser {
    parse(content, property_name = "material") {
        if (property_name !== "material") throw new Error("Implemented for material only");
        return this.parseMaterial();
    }

    parseMaterial() {
        this.cell = this.getCell();
        this.elements = this.getElements();
        this.coordinates = this.getCoordinates();
        this.constraints = this.getConstraints();
        this.units = this.getUnits();
        this.name = this.getName();

        const lattice = Lattice.fromVectors({
            a: this.cell.cell[0],
            b: this.cell.cell[1],
            c: this.cell.cell[2],
            type: this.cell.type,
        });

        const basis = new ConstrainedBasis({
            elements: this.elements,
            coordinates: this.coordinates,
            units: this.units,
            type: this.cell.type,
            cell: lattice.vectorArrays,
            constraints: this.constraints,
        });

        return {
            lattice: lattice.toJSON(),
            basis: basis.toJSON(),
            name: this.name,
            isNonPeriodic: false,
        };
    }

    getCell() {
        throw new Error("Implement in children");
    }

    getElements() {
        throw new Error("Implement in children");
    }

    getCoordinates() {
        throw new Error("Implement in children");
    }

    getConstraints() {
        throw new Error("Implement in children");
    }

    getUnits() {
        throw new Error("Implement in children");
    }

    getName() {
        throw new Error("Implement in children");
    }
}
