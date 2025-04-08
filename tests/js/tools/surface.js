import { Material } from "../../../src/js/material";
import tools from "../../../src/js/tools";
import SiSlab100 from "../../fixtures/si-slab-100.json";
import SiSlab111 from "../../fixtures/si-slab-111.json";
import Si from "../../fixtures/si-standata.json";
import { assertDeepAlmostEqual } from "../utils";

describe("Tools:Surface", () => {
    it("should return slab", () => {
        const material = new Material(Si);
        const slab = tools.surface.generateConfig(material, [1, 0, 0], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab100, slab);
    });

    it("should return slab", () => {
        const material = new Material(Si);
        const slab = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab111, slab);
    });
});
