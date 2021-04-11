import { primitiveCell } from '../../src/cell/primitive_cell';
import { Lattice } from '../../src/lattice/lattice';
import { Na4Cl4 } from '../enums';
import { assertDeepAlmostEqual } from '../utils';

describe('Primitive Cell', () => {
    it('should return primitive lattice', () => {
        const lattice = new Lattice(Na4Cl4.lattice);
        const actualPrimitiveCell = primitiveCell(lattice);
        const expectedPrimitiveCell = [
            [5.691694, 0, 0],
            [0, 5.691694, 0],
            [0, 0, 5.691694],
        ];
        assertDeepAlmostEqual(expectedPrimitiveCell, actualPrimitiveCell);
    });
});
