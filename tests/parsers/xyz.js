import {expect} from "chai";

import {Si} from "../enums";
import {assertDeepAlmostEqual} from "../utils";
import parsers from "../../src/parsers/parsers";

describe('Parsers:XYZ', function () {

    it('should extract basis from XYZ text', function () {
        const text = "Si 0 0 0 \n Si 0.25 0.25 0.25";
        assertDeepAlmostEqual(parsers.xyz.toBasisConfig(text), Si.basis, ["constraints", "cell", "units"]);
    });

});
