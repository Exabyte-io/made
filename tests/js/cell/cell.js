import { Utils } from "@mat3ra/utils";

import { Lattice } from "../../../src/js/lattice/lattice";
import { Silicon } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Cell", () => {
    it("should return scaled cell", () => {
        const lattice = new Lattice(Silicon.lattice);
        const expectedCell = {
            tolerance: 1,
            vector1: [6.697840472868848, 0, 3.867],
            vector2: [2.232613490956283, 6.314784556895033, 3.867],
            vector3: [0, 0, 7.734],
        };
        const actualCell = lattice.Cell.cloneAndScaleByMatrix([
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
        ]);
        assertDeepAlmostEqual(expectedCell, actualCell);
    });
});
