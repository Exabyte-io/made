import { Material } from "../../../src/js/material";
import tools from "../../../src/js/tools";
import SiSlab100 from "../../fixtures/si-slab-100.json";
import SiSlab111 from "../../fixtures/si-slab-111-0.5-vacuum.json";
import SiSlab111NoVacuum from "../../fixtures/si-slab-111-0-vacuum.json";
import Si from "../../fixtures/si-standata.json";
import { assertDeepAlmostEqual } from "../utils";

describe("Tools:Surface", () => {
    it("should return slab", () => {
        const material = new Material(Si);
        const slabConfig = tools.surface.generateConfig(material, [1, 0, 0], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab100, slabConfig);
    });

    it("should return slab", () => {
        const material = new Material(Si);
        const slabConfig = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);

        assertDeepAlmostEqual(SiSlab111NoVacuum, slabConfig);
        const vacuumRatio = 0.5;
        const slabMaterial = new Material(slabConfig);
        const { outOfPlaneAxisIndex } = slabConfig;
        // Add 0.5 vacuum ratio, which is used in MD
        tools.material.scaleOneLatticeVector(
            slabMaterial,
            ["a", "b", "c"][outOfPlaneAxisIndex],
            1 / (1 - vacuumRatio),
        );
        assertDeepAlmostEqual(SiSlab111, slabMaterial.toJSON());
    });
});
