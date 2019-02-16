import {expect} from "chai";

import tools from "../../src/tools";
import {Basis} from "../../src/basis/basis";
import {assertDeepAlmostEqual} from "../utils";
import {Si2Basis, Si2BasisRepeated, H2HFinal, H2HImage, H2HInitial} from "../enums";

describe('Tools:Basis', function () {

    it('should return a repeated basis', function () {
        const b = new Basis(Si2Basis);
        assertDeepAlmostEqual(Si2BasisRepeated, tools.basis.repeat(b, [2, 1, 1]).toJSON());
    });

    it('should return interpolated basises', function () {
        const finalBasis = new Basis(H2HFinal.basis);
        const initialBasis = new Basis(H2HInitial.basis);
        const images = tools.basis.interpolate(initialBasis, finalBasis, 1);
        assertDeepAlmostEqual(H2HImage.basis, images[0].toJSON());
    });

});
