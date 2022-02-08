import { expect } from "chai";

import { Made } from "../../src/made";
import tools from "../../src/tools";
import { Na4Cl4 } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Tools:Material", () => {
    it("should return a material with the 'a' lattice vectors scaled by a factor of 2", () => {
        const material = new Made.Material(Na4Cl4);
        console.log(material);
        tools.material.scaleOneLatticeVector(material, "a", 2.0);
        console.log(material);
        const referenceLatticeA = 11.383388;
        expect(material.lattice.a).to.be.equal(referenceLatticeA);
    });
    it("should return a material with the lattice scaled to the same value twice in a row", () => {
        const material = new Made.Material(Na4Cl4);
        tools.material.scaleLatticeToMakeNonPeriodic(material);
        const updatedLatticeOne = material.lattice;
        tools.material.scaleLatticeToMakeNonPeriodic(material);
        const updatedLatticeTwo = material.lattice;
        assertDeepAlmostEqual(updatedLatticeOne, updatedLatticeTwo);
    });
    it("should return a material with the basis translated to the same center of coordinates twice in a row", () => {
        const material = new Made.Material(Na4Cl4);
        material.lattice = new Made.Lattice(material.lattice);
        tools.material.getBasisConfigTranslatedToCenter(material);
        const updatedBasisOne = material.Basis.coordinates;
        tools.material.getBasisConfigTranslatedToCenter(material);
        const updatedBasisTwo = material.Basis.coordinates;
        assertDeepAlmostEqual(updatedBasisOne, updatedBasisTwo);
    });
    it("should return a material with the basis and lattice scaled to the same values twice in a row", () => {
        const material = new Made.Material(Na4Cl4);
        tools.material.scaleLatticeToMakeNonPeriodic(material);
        const updatedLatticeOne = material.lattice;
        tools.material.getBasisConfigTranslatedToCenter(material);
        const updatedBasisOne = material.Basis.coordinates;
        tools.material.scaleLatticeToMakeNonPeriodic(material);
        const updatedLatticeTwo = material.lattice;
        tools.material.getBasisConfigTranslatedToCenter(material);
        const updatedBasisTwo = material.Basis.coordinates;
        assertDeepAlmostEqual(updatedLatticeOne, updatedLatticeTwo);
        assertDeepAlmostEqual(updatedBasisOne, updatedBasisTwo);
    });
});
