import { expect } from "chai";

import { Material } from "../../src/material";
import { EspressoParser } from "../../src/parsers/espresso/7.2/parser";
import parsers from "../../src/parsers/parsers";
import { NiCub, NiCubPWSCFInput, Si, SiPWSCFInput } from "../enums";

describe("Parsers:Espresso", () => {
    it("should return textual representation of a material according to QE pw.x input format", () => {
        const material = new Material(Si);
        expect(parsers.espresso.toEspressoFormat(material)).to.be.equal(SiPWSCFInput);
    });

    it("should return inetrmediate format for Ni ", () => {
        const parser = new EspressoParser(NiCubPWSCFInput);
        expect(parser.getIntermediateFormat()).to.be.equal(NiCub);
    });
});
