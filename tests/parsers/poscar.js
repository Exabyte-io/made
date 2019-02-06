import {expect} from "chai";

import {Material} from "../../src/material";
import {Na4Cl4, Na4Cl4Poscar} from "../enums";

describe('POSCAR', function () {
    it('should return poscar with sorted basis items', function () {
        const material = new Material(Na4Cl4);
        expect(`${material.getAsPOSCAR()}\n`).to.be.equal(Na4Cl4Poscar);
    });
});
