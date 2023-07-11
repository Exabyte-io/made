import { expect } from "chai";
import { mix } from "mixwith";

import { BaseParser } from "../../src/parsers/init";
import { FortranParserMixin } from "../../src/parsers/utils/fortran";
import { FortranFile1, FortranFile1JSON } from "../enums";

// describe("Parsers:Fortran", () => {
//     class TestParser extends mix(BaseParser).with(FortranParserMixin) {} // Test class
//     it("should return intermediate format of parsed input file", () => {
//         const parser = new TestParser({});
//         const data = parser.parse(FortranFile1);
//         expect(data).to.be.deep.equal(FortranFile1JSON);
//     });
// });
