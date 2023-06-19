import { expect } from "chai";

import { parseFortranFile } from "../../src/parsers/fortran/fortran";
import { FortranFile1, FortranFile1JSON, FortranFileInvalid } from "../enums";

describe("parseFortranFile", () => {
    it("should return json for FortranFile1", () => {
        expect(parseFortranFile(FortranFile1)).to.be.deep.equal(FortranFile1JSON);
    });
    it("should throw error for non-Fortran file", () => {
        expect(() => parseFortranFile(FortranFileInvalid)).to.throw("Incorrect fortran file");
    });
});
