import { expect } from "chai";

import { parseFortranFile } from "../../src/parsers/fortran_parser";
import { FortranFile1, FortranFile1JSON } from "../enums";

describe("parseFortranFile", () => {
    it("should return json for FortranFile1", () => {
        expect(parseFortranFile(FortranFile1)).to.be.deep.equal(FortranFile1JSON);
    });
});
