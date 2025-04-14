import { Utils } from "@mat3ra/utils";

import { Material } from "../../../src/js/material";
import tools from "../../../src/js/tools";
import { Silicon, SiSlab100, SiSlab111, SiSlab111NoVacuum } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Tools:Surface", () => {
    it("should return slab (100)", () => {
        const material = new Material(Silicon);
        const slabConfig = tools.surface.generateConfig(material, [1, 0, 0], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab100, slabConfig);
    });

    it("should return slab (111)", () => {
        const material = new Material(Silicon);
        const slabConfig = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab111NoVacuum, slabConfig);
    });

    it("should return slab (111) with vacuum", () => {
        const material = new Material(Silicon);
        const slabConfig = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);
        const vacuumRatio = 0.5;
        const slabMaterial = new Material(slabConfig);
        const { outOfPlaneAxisIndex } = slabConfig;
        // Add 0.5 vacuum ratio, which is used in MD
        tools.material.scaleOneLatticeVector(
            slabMaterial,
            ["a", "b", "c"][outOfPlaneAxisIndex],
            1 / (1 - vacuumRatio),
        );
        const expectedSlabMaterial = new Material(SiSlab111);
        assertDeepAlmostEqual(expectedSlabMaterial.toJSON(), slabMaterial.toJSON());
    });
});
