import { Lattice } from "../../src/lattice/lattice";
import { Material } from "../../src/material";
import tools from "../../src/tools";
import { C2H4, C2H4Translated, Na4Cl4 } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Tools:Material", () => {
    it("should return a scaled lattice vector", () => {
        const newMaterial = new Material(Na4Cl4);
        const latticeVectorKey = "a";
        const factor = 2.0;
        const expectedLatticeVectorArray = [
            newMaterial.lattice.a * factor,
            newMaterial.lattice.b,
            newMaterial.lattice.c,
        ];
        tools.material.scaleOneLatticeVector(newMaterial, latticeVectorKey, factor);
        const actualLatticeVectorArray = [
            newMaterial.lattice.a,
            newMaterial.lattice.b,
            newMaterial.lattice.c,
        ];
        assertDeepAlmostEqual(expectedLatticeVectorArray, actualLatticeVectorArray);
    });
    it("should return a scaled non-periodic lattice", () => {
        const newMaterial = new Material(C2H4);
        const expectedLatticeValue = 3.162;
        tools.material.scaleLatticeToMakeNonPeriodic(newMaterial);
        assertDeepAlmostEqual(expectedLatticeValue, newMaterial.lattice.a);
        assertDeepAlmostEqual(expectedLatticeValue, newMaterial.lattice.b);
        assertDeepAlmostEqual(expectedLatticeValue, newMaterial.lattice.c);
    });
    it("should return a material basis translated to the center of the material lattice", () => {
        const newMaterial = new Material(C2H4);
        newMaterial.lattice = new Lattice(newMaterial.lattice);
        tools.material.getBasisConfigTranslatedToCenter(newMaterial);
        assertDeepAlmostEqual(newMaterial.basis.coordinates, C2H4Translated.basis.coordinates);
    });
});
