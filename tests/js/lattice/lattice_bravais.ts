import { Utils } from "@mat3ra/utils";
import { expect } from "chai";

import { Lattice } from "../../../src/js/lattice/lattice";
import { Na4Cl4 } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Lattice Bravais", () => {
    it("should return lattice bravais from vectors", () => {
        const lattice = Lattice.fromVectors({
            ...Na4Cl4.lattice.vectors,
            type: "CUB",
        });
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice, ["type", "vectors"]);
    });

    it("should return a list of editable keys", () => {
        const lattice = new Lattice.fromConfig({
            ...Na4Cl4.lattice,
            type: "CUB",
        });
        expect(lattice.editables).to.be.deep.equal({ a: true });
    });

    it("lattice created from vectors should be equal to the original lattice", () => {
        const lattice = Lattice.fromVectors({
            ...Na4Cl4.lattice.vectors,
            type: "CUB",
        });
        const originalLattice = new Lattice({
            ...Na4Cl4.lattice,
            type: "CUB",
        });

        assertDeepAlmostEqual(lattice.toJSON(), originalLattice.toJSON(), ["type"]);
    });
});
