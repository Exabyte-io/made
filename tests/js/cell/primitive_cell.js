import { Utils } from "@mat3ra/utils";

import { primitiveCell } from "../../../src/js/cell/primitive_cell";
import { Lattice } from "../../../src/js/lattice/lattice";
import { Na4Cl4 } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Primitive Cell", () => {
    it("should return primitive lattice", () => {
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
