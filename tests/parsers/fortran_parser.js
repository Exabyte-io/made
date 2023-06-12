import { expect } from "chai";

import { parseFortranFile } from "../../src/parsers/fortran_parser";
import { FortranFile1, FortranFile1JSON, FortranFileInvalid, FortranFileNoCards } from "../enums";

describe("parseFortranFile", () => {
    it("should return json for FortranFile1", () => {
        expect(parseFortranFile(FortranFile1)).to.be.deep.equal(FortranFile1JSON);
    });
    it("should throw error for non-Fortran file", () => {
        expect(() => parseFortranFile(FortranFileInvalid)).to.throw("Incorrect fortran file");
    });
    it("should throw error for Fortran file without cards", () => {
        expect(() => parseFortranFile(FortranFileNoCards)).to.throw("No cards found");
    });
});
