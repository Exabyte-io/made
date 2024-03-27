import { expect } from "chai";

import { Basis } from "../../../src/js/basis/basis";
import { Made } from "../../../src/js/made";
import { Material } from "../../../src/js/material";
import { Si, SiSupercell } from "../enums";

describe("Tools:Supercell", () => {
    it("should generate supercell", () => {
        const material = new Material(Si);
        const supercell = Made.tools.supercell.generateConfig(material, [
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
        ]);
        expect(supercell.lattice).deep.equal(SiSupercell.lattice, "lattices are not equal");

        const basis1 = new Basis(SiSupercell.basis);
        const basis2 = new Basis(supercell.basis);

        expect(basis1.isEqualTo(basis2)).to.be.equal(true);
    });
});
