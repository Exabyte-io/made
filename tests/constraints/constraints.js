import {expect} from "chai";

import {atomicConstraints} from "../enums";
import parsers from "../../src/parsers/parsers";
import {AtomicConstraints} from "../../src/constraints/constraints";

describe('AtomicConstraints', function () {

    it('should return constraints as string', function () {
        const constraints = new AtomicConstraints(atomicConstraints);
        expect(constraints.getAsStringByIndex(0)).to.be.equal("1 1 0");
    });

    it('should return constraints as string with given map function', function () {
        const constraints = new AtomicConstraints(atomicConstraints);
        expect(constraints.getAsStringByIndex(0, parsers.poscar.atomicConstraintsCharFromBool)).to.be.equal("T T F");
    });

});
