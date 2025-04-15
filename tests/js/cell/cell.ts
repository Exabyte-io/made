import { Utils } from "@mat3ra/utils";

import { Lattice } from "../../../src/js/lattice/lattice";
import { Silicon } from "../fixtures";
import { LatticeSchema } from "@mat3ra/esse/dist/js/types";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Cell", () => {
    it("should return scaled cell", () => {
        const lattice = new Lattice(Silicon.lattice as LatticeSchema);
        const expectedCell = {
            a: [6.697840472868848, 0, 3.867],
            b: [2.232613490956283, 6.314784556895033, 3.867],
            c: [0, 0, 7.734],
        };
        const actualCell = lattice.Cell.cloneAndScaleByMatrix([
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
        ]);
        assertDeepAlmostEqual(expectedCell, actualCell);
    });
});
