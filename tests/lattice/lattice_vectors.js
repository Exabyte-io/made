import {expect} from "chai";

import {Na4Cl4} from "../enums";
import {assertDeepAlmostEqual} from "../utils";
import {LatticeVectors} from "../../src/lattice/lattice_vectors";

describe('Lattice Vectors', function () {

    it('should return lattice from bravais', function () {
        const lattice = LatticeVectors.fromBravais(Na4Cl4.lattice);
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice.vectors, ["type", "vectors"]);
    });

});
