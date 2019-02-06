import {expect} from "chai";

import {Na4Cl4} from "../enums";
import {Basis} from "../../src/basis/basis";

describe('Basis', function () {

    it('should return proper unitCellFormula', function () {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.unitCellFormula).to.be.equal('Na4Cl4');
    });

    it('should return proper formula', function () {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.formula).to.be.equal('NaCl');
    });

});
