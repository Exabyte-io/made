import { expect } from "chai";

import nativeFormatParsers from "../../src/parsers/native_format_parsers";
import {
    BNHex,
    BNHexIbravPWSCFInput,
    BNHexPWSCFInput,
    Graphene,
    GraphenePoscar,
    GraphenePWSCFInput,
    NiCub,
    NiCubAPWSCFInput,
    NiCubCPPWSCFInput,
    NiCubPWSCFInput,
    NiHex,
    NiHexPoscar,
    Sb2S3Orc,
    Sb2S3OrcAPWSCFInput,
    Sb2S3OrcPWSCFInput,
} from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers.NativeFormat", () => {
    it("should return a material config for graphene from a json", () => {
        const json = JSON.stringify(Graphene);
        expect(nativeFormatParsers.convertFromNativeFormat(json)).to.deep.equal(Graphene);
    });

    it("should return a material config for graphene from a poscar", () => {
        const poscar = GraphenePoscar;
        const config = nativeFormatParsers.getMaterialConfig(poscar);
        expect(config).to.be.deep.almost.equal(Graphene);
    });

    it("should return a material config for Ni hex from a poscar", () => {
        const poscar = NiHexPoscar;
        const config = nativeFormatParsers.getMaterialConfig(poscar);
        assertDeepAlmostEqual(config, NiHex, ["lattice"]);
        assertDeepAlmostEqual(config.lattice, NiHex.lattice, ["type"]); // to omit "lattice.type" property
    });

    it("should return a material config for graphene from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(GraphenePWSCFInput);
        assertDeepAlmostEqual(config, Graphene, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(Graphene.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for Ni CUB with specified CELL_PARAMETERS from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(NiCubCPPWSCFInput);
        assertDeepAlmostEqual(config, NiCub, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(NiCub.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for Ni CUB with specified parameter celldm(1) from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(NiCubPWSCFInput);
        assertDeepAlmostEqual(config, NiCub, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(NiCub.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for Ni CUB with specified parameter A from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(NiCubAPWSCFInput);
        assertDeepAlmostEqual(config, NiCub, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(NiCub.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for Sb2S3 ORC from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(Sb2S3OrcPWSCFInput);
        assertDeepAlmostEqual(config, Sb2S3Orc, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(Sb2S3Orc.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for Sb2S3 ORC with specified parameter A from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(Sb2S3OrcAPWSCFInput);
        assertDeepAlmostEqual(config, Sb2S3Orc, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(Sb2S3Orc.name.toLowerCase()); // to compare case insensitively
    });

    it("should return a material config for BN HEX from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(BNHexPWSCFInput);
        assertDeepAlmostEqual(config, BNHex, ["name", "lattice"]); // title is omitted in input file
        assertDeepAlmostEqual(config.lattice, BNHex.lattice, ["type"]); // lattice type isn't detected currently
    });

    it("should return a material config for BN HEX with specified ibrav from a QE input file", () => {
        const config = nativeFormatParsers.getMaterialConfig(BNHexIbravPWSCFInput);
        assertDeepAlmostEqual(config, BNHex, ["name"]); // title is omitted in input file
    });

    it("should throw an error for unknown format", () => {
        const text = "A\n snippet from an unknown format";
        expect(() => nativeFormatParsers.getMaterialConfig(text)).to.throw("Unknown format");
    });
});
