import { Utils } from "@mat3ra/utils";

import { Basis } from "../../../src/js/basis/basis";
import tools from "../../../src/js/tools";
import { H2HFinal, H2HImage, H2HInitial, Si2Basis, Si2BasisRepeated } from "../fixtures";

const { assertDeepAlmostEqual } = Utils.assertion;

describe("Tools:Basis", () => {
    it("should return a repeated basis", () => {
        const basis = new Basis(Si2Basis);
        const repeatedBasis = tools.basis.repeat(basis, [2, 1, 1]);
        assertDeepAlmostEqual(Si2BasisRepeated, repeatedBasis.toJSON());
    });

    it("should return interpolated basises", () => {
        const finalBasis = new Basis(H2HFinal.basis);
        const initialBasis = new Basis(H2HInitial.basis);
        const images = tools.basis.interpolate(initialBasis, finalBasis, 1);
        assertDeepAlmostEqual(H2HImage.basis, images[0].toJSON());
    });
});
