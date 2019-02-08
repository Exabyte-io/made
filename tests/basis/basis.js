import {expect} from "chai";

import {Basis} from "../../src/basis/basis";
import {FeLiSiBasis, LiFeSiBasis, Na4Cl4} from "../enums";

describe('Basis', function () {

    it('should return same basis', function () {
        const basis1 = new Basis(FeLiSiBasis);
        const basis2 = new Basis(LiFeSiBasis);
        expect(basis1.isEqualTo(basis2)).to.be.ok;
    });

    it('should return proper unitCellFormula', function () {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.unitCellFormula).to.be.equal('Na4Cl4');
    });

    it('should return proper formula', function () {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.formula).to.be.equal('NaCl');
    });

});
