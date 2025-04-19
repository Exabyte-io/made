import { Utils } from "@mat3ra/utils";

import { Material } from "../../../src/js/material";
import tools from "../../../src/js/tools";
import { Silicon, SiSlab100, SiSlab111, SiSlab111Gamma120, SiSlab111NoVacuum } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

const generateSlabWithVacuum = (slabConfig, vacuumRatio) => {
    const slabMaterial = new Material(slabConfig);
    const { outOfPlaneAxisIndex } = slabConfig;

    tools.material.scaleOneLatticeVector(
        slabMaterial,
        ["a", "b", "c"][outOfPlaneAxisIndex],
        1 / (1 - vacuumRatio),
    );

    return slabMaterial;
};

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
        const slabMaterial = generateSlabWithVacuum(slabConfig, vacuumRatio);
        const expectedSlabMaterial = new Material(SiSlab111);
        const expectedMaterialJSON = expectedSlabMaterial.toJSON();
        const slabMaterialJSON = slabMaterial.toJSON();
        assertDeepAlmostEqual(expectedMaterialJSON, slabMaterialJSON);
    });

    it("should return slab (111) with vacuum (gamma~=120) for gamma = 59.999", () => {
        const adjustedSilicon = {
            ...Silicon,
            lattice: {
                ...Silicon.lattice,
                gamma: 59.999,
            },
        };
        const material = new Material(adjustedSilicon);
        const slabConfig = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);
        const vacuumRatio = 0.5;
        const slabMaterial = generateSlabWithVacuum(slabConfig, vacuumRatio);
        // With original gamma being below 60, the surface generated differently, with a different gamma
        const expectedSlabMaterial = new Material(SiSlab111Gamma120);
        const expectedMaterialJSON = expectedSlabMaterial.toJSON();
        const slabMaterialJSON = slabMaterial.toJSON();
        assertDeepAlmostEqual(expectedMaterialJSON, slabMaterialJSON);
    });

    it("should return slab (111) with vacuum for gamma = 60.001", () => {
        const adjustedSilicon = {
            ...Silicon,
            lattice: {
                ...Silicon.lattice,
                gamma: 60.001,
            },
        };
        const material = new Material(adjustedSilicon);
        const slabConfig = tools.surface.generateConfig(material, [1, 1, 1], 3, 1, 1);
        const vacuumRatio = 0.5;
        const slabMaterial = generateSlabWithVacuum(slabConfig, vacuumRatio);
        const expectedSlabMaterial = new Material(SiSlab111);
        const expectedMaterialJSON = expectedSlabMaterial.toJSON();
        const slabMaterialJSON = slabMaterial.toJSON();
        assertDeepAlmostEqual(expectedMaterialJSON, slabMaterialJSON);
    });
});
