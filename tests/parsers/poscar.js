import {expect} from "chai";

import {Material} from "../../src/material";
import {Na4Cl4, Na4Cl4Poscar, Zr12H24, Zr12H24Poscar} from "../enums";

describe('Parsers.POSCAR', function () {
    it('should return poscar with sorted basis items', function () {
        const material = new Material(Na4Cl4);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Na4Cl4Poscar);
    });

    it('should return poscar with "H Zr" in elements line', function () {
        const material = new Material(Zr12H24);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Zr12H24Poscar);
    });

});
