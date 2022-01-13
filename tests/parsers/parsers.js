import { xyzToPoscar } from "../../src/parsers/parsers";
import { CH4, CH4POSCAR } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Parsers:XYZ", () => {
    it("should return the xyz file content in poscar file format", () => {
        assertDeepAlmostEqual(xyzToPoscar(CH4), CH4POSCAR);
    });
});
