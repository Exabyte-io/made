import {Basis} from "./basis";
import {expect} from "chai";

/*
 * @summary Tests boilerplate. Can be run as below:
 * ./node_modules/.bin/mocha --compilers js:node_modules/babel-core/register.js imports/made/basis/basis.tests.js
 */
describe('Basis', function () {
    it('isEqualTo', function () {
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
});
