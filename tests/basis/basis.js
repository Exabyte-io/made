import {expect} from "chai";

import {Na4Cl4} from "../enums";
import {Basis} from "../../src/basis/basis";

describe('Basis', function () {

    it('should return same basis', function () {
        const basis1 = new Basis({
            elements: [
                {
                    id: 1,
                    value: 'Fe'
                },
                {
                    id: 2,
                    value: 'Li'
                },
                {
                    id: 3,
                    value: 'Si'
                }
            ],
            coordinates: [
                {
                    id: 1,
                    value: [0.5, 0.5, 0.5]
                },
                {
                    id: 2,
                    value: [0.3, 0.3, 0.3]
                },
                {
                    id: 3,
                    value: [0, 0, 0]
                }
            ]
        });
        const basis2 = new Basis({
            elements: [
                {
                    id: 2,
                    value: 'Fe'
                },
                {
                    id: 1,
                    value: 'Li'
                },
                {
                    id: 3,
                    value: 'Si'
                }
            ],
            coordinates: [
                {
                    id: 2,
                    value: [0.5, 0.5, 0.5]
                },
                {
                    id: 1,
                    value: [0.3, 0.3, 0.3]
                },
                {
                    id: 3,
                    value: [0, 0, 0]
                }
            ]
        });
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
