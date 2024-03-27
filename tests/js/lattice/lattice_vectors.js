import { LatticeVectors } from "../../../src/js/lattice/lattice_vectors";
import { Na4Cl4 } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Lattice Vectors", () => {
    it("should return lattice from bravais", () => {
        const lattice = LatticeVectors.fromBravais(Na4Cl4.lattice);
        assertDeepAlmostEqual(lattice.toJSON(), Na4Cl4.lattice.vectors, ["type", "vectors"]);
    });
});
