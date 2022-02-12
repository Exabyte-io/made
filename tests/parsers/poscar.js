import { expect } from "chai";

import { Material } from "../../src/material";
import { convertFromOtherFormat, getAtomsCount } from "../../src/parsers/poscar";
import {
    CH4,
    CH4POSCAR,
    H2O,
    Na4Cl4,
    Na4Cl4Poscar,
    Zr1H23Zr1H1,
    Zr1H23Zr1H1Poscar,
} from "../enums";
import { assertDeepAlmostEqual } from "../utils";

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
        expect(getAtomsCount(H2O)).to.be.equal(3);
    });
    it("should return the xyz file content in poscar file format", () => {
        assertDeepAlmostEqual(convertFromOtherFormat(CH4, "xyz"), CH4POSCAR);
    });
});
