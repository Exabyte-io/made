import {expect} from "chai";

import {Si, SiSlab} from "../enums";
import tools from "../../src/tools";
import {Material} from "../../src/material";
import {assertDeepAlmostEqual} from "../utils";

describe('Tools:Surface', function () {

    it('should return slab', function () {
        const material = new Material(Si);
        const slab = tools.surface.generateConfig(material, [1, 0, 0], 3, 1, 1);
        assertDeepAlmostEqual(SiSlab, slab);
    });

});
