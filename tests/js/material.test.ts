import { expect } from "chai";

import { Material } from "../../src/js/material";
import { Na4Cl4 } from "./fixtures";

describe("Material", () => {
    it("should return unique elements", () => {
        const material = new Material(Na4Cl4);
        expect(material.uniqueElements).to.have.same.members(["Na", "Cl"]);
    });
});
