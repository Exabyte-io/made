import {expect} from "chai";

import {Na4Cl4, Si} from "../enums";
import {assertDeepAlmostEqual} from "../utils";
import {Lattice} from "../../src/lattice/lattice";

describe('Lattice', function () {

    it('should return lattice cell volume', function () {
        const lattice = new Lattice(Si.lattice);
        expect(lattice.volume).to.be.almost.equal(125);
    });

    it('should return lattice from vectors', function () {
        const lattice = Lattice.fromVectors(Na4Cl4.lattice.vectors);
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice, ["type"])
    });

    it('should return lattice type', function () {
        const lattice = new Lattice(Si.lattice);
        expect(lattice.typeExtended).to.be.equal("TRI_1b");
    });

    /**
     * hash
     */

    it('should return lattice hash string', function () {
        const lattice = new Lattice(Na4Cl4.lattice);
        expect(lattice.getHashString()).to.be.equal('5.692;5.692;5.692;90;90;90;');
    });
});
