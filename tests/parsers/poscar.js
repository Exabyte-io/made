import {expect} from "chai";

import {Material} from "../../src/material";
import {atomsCount} from "../../src/parsers/poscar";
import {Na4Cl4, Na4Cl4Poscar, Zr1H23Zr1H1, Zr1H23Zr1H1Poscar, H2O} from "../enums";

describe('Parsers.POSCAR', function () {
    it('should return a valid poscar', function () {
        const material = new Material(Na4Cl4);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Na4Cl4Poscar);
    });

    it('should return poscar elements line according to id in basis, duplicate entries separate', function () {
        const material = new Material(Zr1H23Zr1H1);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Zr1H23Zr1H1Poscar);
    });

    it('should return the number of atoms for a molecule in a poscar file', function() {
        expect(atomsCount(H2O)).to.be.equal(3);
    });

});
