import {expect} from "chai";

import {Made} from "../../src/made";
import {Si, SiSupercell} from "../enums";
import {Basis} from "../../src/basis/basis";
import {Material} from "../../src/material";

describe('Tools:Supercell', function () {

    it('should generate supercell', function () {
        const material = new Material(Si);

        const supercell = Made.tools.supercell.generateConfig(material, [[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        expect(supercell.lattice).deep.equal(SiSupercell.lattice, 'lattices are not equal');

        const basis1 = new Basis(SiSupercell.basis);
        const basis2 = new Basis(supercell.basis);

        expect(basis1.isEqualTo(basis2)).to.be.equal(true);
    });
});
