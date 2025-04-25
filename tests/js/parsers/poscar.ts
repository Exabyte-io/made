import { expect } from "chai";

import { AtomicConstraints } from "../../../src/js/constraints/constraints";
import { Material } from "../../../src/js/material";
import parsers from "../../../src/js/parsers/parsers";
import { atomsCount } from "../../../src/js/parsers/poscar";
import {
    atomicConstraints,
    H2OPoscar,
    Na4Cl4,
    Na4Cl4Poscar,
    Zr1H23Zr1H1,
    Zr1H23Zr1H1Poscar,
} from "../fixtures";

describe("Parsers.POSCAR", () => {
    it("should return a valid poscar", () => {
        const material = new Material(Na4Cl4);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Na4Cl4Poscar);
    });

    it("should return poscar elements line according to id in basis, duplicate entries separate", () => {
        const material = new Material(Zr1H23Zr1H1);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Zr1H23Zr1H1Poscar);
    });

    it("should return the number of atoms for a molecule in a poscar file", () => {
        expect(atomsCount(H2OPoscar)).to.be.equal(3);
    });

    it("should return constraints as string with given map function", () => {
        const constraints = AtomicConstraints.fromObjects(atomicConstraints);
        expect(
            constraints.getAsStringByIndex(0, parsers.poscar.atomicConstraintsCharFromBool),
        ).to.be.equal("T T F");
    });
});
