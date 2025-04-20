import "../setup";

import { LatticeSchema, LatticeVectorsSchema } from "@mat3ra/esse/dist/js/types";
import { Utils } from "@mat3ra/utils";
import { expect } from "chai";

import { Lattice } from "../../../src/js/lattice/lattice";
import { Na4Cl4, Silicon } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Lattice", () => {
    it("should create a lattice", () => {
        const lattice = Lattice.fromConfig(Silicon.lattice);
        expect(lattice).to.be.instanceOf(Lattice);
        expect(lattice.a).to.be.equal(Silicon.lattice.a);
    });
    it("should return lattice type", () => {
        const lattice = new Lattice(Silicon.lattice as LatticeSchema);
        expect(lattice.type).to.be.equal("FCC");
    });

    it("should return lattice cell volume", () => {
        const lattice = new Lattice(Silicon.lattice as LatticeSchema);
        expect(lattice.volume).to.be.almost.equal(40.889096881496656);
    });

    it("should return lattice from vectors", () => {
        const lattice = Lattice.fromVectors(Na4Cl4.lattice.vectors as LatticeVectorsSchema);
        assertDeepAlmostEqual(Na4Cl4.lattice, lattice.toJSON(), ["type"]);
    });

    it("should return lattice type", () => {
        const lattice = new Lattice(Silicon.lattice as LatticeSchema);
        expect(lattice.typeExtended).to.be.equal("FCC");
    });

    it("should return lattice hash string", () => {
        const lattice = new Lattice(Na4Cl4.lattice as LatticeSchema);
        expect(lattice.getHashString()).to.be.equal("5.692;5.692;5.692;90;90;90;");
    });
    // TODO: add more tests
});
