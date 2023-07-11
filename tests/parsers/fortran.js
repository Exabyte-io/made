import { expect } from "chai";

import { FortranParser } from "../../src/parsers/utils/fortran";
import { FortranFile1, FortranFile1JSON } from "../enums";

describe("Parsers:Fortran", () => {
    it("should return intermediate format of parsed input file", () => {
        const parser = new FortranParser();
        const data = parser.parse(FortranFile1);
        expect(data).to.be.equal(FortranFile1JSON);
    });
});
