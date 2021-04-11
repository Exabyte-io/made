import { expect } from 'chai';

import { AtomicConstraints } from '../../src/constraints/constraints';
import parsers from '../../src/parsers/parsers';
import { atomicConstraints } from '../enums';
import { assertDeepAlmostEqual } from '../utils';

describe('AtomicConstraints', () => {
    it('should return constraints as string', () => {
        const constraints = new AtomicConstraints(atomicConstraints);
        expect(constraints.getAsStringByIndex(0)).to.be.equal('1 1 0');
    });

    it('should return constraints as string with given map function', () => {
        const constraints = new AtomicConstraints(atomicConstraints);
        expect(
            constraints.getAsStringByIndex(0, parsers.poscar.atomicConstraintsCharFromBool),
        ).to.be.equal('T T F');
    });

    it('should return [true, true, true] as constraints for the line without any constraints', () => {
        const text = 'Si 0 0 0 0 1 1\n Si 0.25 0.25 0.25';
        const basis = parsers.xyz.toBasisConfig(text);
        assertDeepAlmostEqual(basis.constraints[0].value, [false, true, true]);
        assertDeepAlmostEqual(basis.constraints[1].value, [true, true, true]);
    });
});
