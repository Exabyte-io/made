import {expect} from "chai";

import {Na4Cl4} from "../enums";
import {assertDeepAlmostEqual} from "../utils";
import {LatticeBravais} from "../../src/lattice/lattice_bravais";

describe('Lattice Bravais', function () {

    it('should return lattice bravais from vectors', function () {
        const lattice = LatticeBravais.fromVectors(Na4Cl4.lattice.vectors);
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice, ["type", "vectors"]);
    });

    it('should return a list of editable keys', function () {
        const lattice = new LatticeBravais(Na4Cl4.lattice);
        expect(lattice.editables).to.be.deep.equal({a: true});
    });

});
