import parsers from '../../src/parsers/parsers';
import { Si } from '../enums';
import { assertDeepAlmostEqual } from '../utils';

describe('Parsers:XYZ', () => {
    it('should extract basis from XYZ text', () => {
        const text = 'Si 0 0 0 \n Si 0.25 0.25 0.25';
        assertDeepAlmostEqual(parsers.xyz.toBasisConfig(text), Si.basis, [
            'constraints',
            'cell',
            'units',
        ]);
    });
});
