import { Utils } from "@mat3ra/utils";
import { expect } from "chai";

import { LatticeBravais } from "../../../src/js/lattice/lattice_bravais";
import { Na4Cl4 } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Lattice Bravais", () => {
    it("should return lattice bravais from vectors", () => {
        const lattice = LatticeBravais.fromVectors(Na4Cl4.lattice.vectors);
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice, ["type", "vectors"]);
    });

    it("should return a list of editable keys", () => {
        const lattice = new LatticeBravais(Na4Cl4.lattice);
        expect(lattice.editables).to.be.deep.equal({ a: true });
    });

    it("lattice created from vectors should be equal to the original lattice", () => {
        const lattice = LatticeBravais.fromVectors(Na4Cl4.lattice.vectors);
        const originalLattice = new LatticeBravais(Na4Cl4.lattice);

        assertDeepAlmostEqual(lattice.toJSON(), originalLattice.toJSON(), ["type"]);
    });
});
