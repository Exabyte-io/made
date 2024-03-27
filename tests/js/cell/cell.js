import { Lattice } from "../../../src/js/lattice/lattice";
import { Si } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Cell", () => {
    it("should return scaled cell", () => {
        const lattice = new Lattice(Si.lattice);
        const expectedCell = {
            tolerance: 1,
            vector1: [10, 0, 0],
            vector2: [0, 10, 0],
            vector3: [0, 0, 10],
        };
        const actualCell = lattice.Cell.cloneAndScaleByMatrix([
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
        ]);
        assertDeepAlmostEqual(expectedCell, actualCell);
    });
});
