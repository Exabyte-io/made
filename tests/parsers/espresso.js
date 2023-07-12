import { expect } from "chai";

import { Material } from "../../src/material";
import { ESPRESSOMaterialParser } from "../../src/parsers/espresso/parser";
import parsers from "../../src/parsers/parsers";
import {
    BNHex,
    BNHexIbravPWSCFInput,
    BNHexPWSCF,
    NiCub,
    NiCubIbravAPWSCFInput,
    Si,
    SiPWSCFInput,
} from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers:Espresso", () => {
    it("should return textual representation of a material according to QE pw.x input format", () => {
        const material = new Material(Si);
        expect(parsers.espresso.toEspressoFormat(material)).to.be.equal(SiPWSCFInput);
    });

    it("should return a material config from QE input file for BN Hex with specified ibrav and celldm parameter", () => {
        const parser = new ESPRESSOMaterialParser();
        const materialConfig = parser.parse(BNHexIbravPWSCFInput, "material");
        assertDeepAlmostEqual(materialConfig, BNHex, ["name"]);
    });

    it("should return a material config from QE input file for BN Hex with cell parameters given", () => {
        const parser = new ESPRESSOMaterialParser();
        const materialConfig = parser.parse(BNHexPWSCF, "material");
        assertDeepAlmostEqual(materialConfig, BNHex, ["name", "lattice"]); // lattice.type is not detected, defaults to TRI, skipping it in tests
        assertDeepAlmostEqual(materialConfig.lattice, BNHex.lattice, ["type"]);
    });

    it("should return a material config from QE input file for Ni Cub with specified ibrav and A parameter", () => {
        const parser = new ESPRESSOMaterialParser();
        const materialConfig = parser.parse(NiCubIbravAPWSCFInput, "material");
        assertDeepAlmostEqual(materialConfig, NiCub, ["name"]);
    });
});
