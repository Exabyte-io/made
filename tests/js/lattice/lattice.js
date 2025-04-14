import { Utils } from "@mat3ra/utils";
import { expect } from "chai";

import { Lattice } from "../../../src/js/lattice/lattice";
import { Na4Cl4, Silicon } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Lattice", () => {
    it("should return lattice cell volume", () => {
        const lattice = new Lattice(Silicon.lattice);
        expect(lattice.volume).to.be.almost.equal(40.889096881496656);
    });

    it("should return lattice from vectors", () => {
        const lattice = Lattice.fromVectors(Na4Cl4.lattice.vectors);
        assertDeepAlmostEqual(Na4Cl4.lattice, lattice.toJSON(), ["type"]);
    });

    it("should return lattice type", () => {
        const lattice = new Lattice(Silicon.lattice);
        expect(lattice.typeExtended).to.be.equal("FCC");
    });

    it("should return lattice hash string", () => {
        const lattice = new Lattice(Na4Cl4.lattice);
        expect(lattice.getHashString()).to.be.equal("5.692;5.692;5.692;90;90;90;");
    });
});
