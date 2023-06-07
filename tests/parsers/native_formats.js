import { expect } from "chai";

import nativeFormatParsers from "../../src/parsers/native_format_parsers";
import {
    Graphene,
    GraphenePoscar,
    GraphenePWSCFInput,
    NiCub,
    NiCubPWSCFInput,
    NiHex,
    NiHexPoscar,
} from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers.NativeFormat", () => {
    it("should return a material config for graphene from a json", () => {
        const json = JSON.stringify(Graphene);
        expect(nativeFormatParsers.convertFromNativeFormat(json)).to.deep.equal(Graphene);
    });

    it("should return a material config for graphene from a poscar", () => {
        const poscar = GraphenePoscar;
        const config = nativeFormatParsers.convertFromNativeFormat(poscar);
        expect(config).to.be.deep.almost.equal(Graphene);
    });

    it("should return a material config for Ni hex from a poscar", () => {
        const poscar = NiHexPoscar;
        const config = nativeFormatParsers.convertFromNativeFormat(poscar);
        assertDeepAlmostEqual(config, NiHex, ["lattice"]);
        assertDeepAlmostEqual(config.lattice, NiHex.lattice, ["type"]); // to omit "lattice.type" property
    });

    it("should return a material config for graphene from a QE input file", () => {
        const qein = GraphenePWSCFInput;
        const config = nativeFormatParsers.convertFromNativeFormat(qein);
        assertDeepAlmostEqual(config, Graphene, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(Graphene.name.toLowerCase());
    });

    it("should return a material config for Ni cub from a QE input file", () => {
        const qein = NiCubPWSCFInput;
        const config = nativeFormatParsers.convertFromNativeFormat(qein);
        console.log("RESULT:", config);
        console.log("EXPECTED:", NiCub);
        assertDeepAlmostEqual(config, NiCub, ["name"]);
        expect(config.name.toLowerCase()).to.be.equal(NiCub.name.toLowerCase());
    });

    it("should throw an error for unknown format", () => {
        const text = "A\n snippet from an unknown format";
        expect(() => nativeFormatParsers.convertFromNativeFormat(text)).to.throw("Unknown format");
    });
});
