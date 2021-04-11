import { Material } from '../../src/material';
import tools from '../../src/tools';
import { Si, SiSlab } from '../enums';
import { assertDeepAlmostEqual } from '../utils';

describe('Tools:Surface', () => {
    it('should return slab', () => {
        const material = new Material(Si);
        const slab = tools.surface.generateConfig(material, [1, 0, 0], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab, slab);
    });
});
