import { Basis } from "../../../src/js/basis/basis";
import tools from "../../../src/js/tools";
import { H2HFinal, H2HImage, H2HInitial, Si2Basis, Si2BasisRepeated } from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Tools:Basis", () => {
    it("should return a repeated basis", () => {
        const b = new Basis(Si2Basis);
        assertDeepAlmostEqual(Si2BasisRepeated, tools.basis.repeat(b, [2, 1, 1]).toJSON());
    });

    it("should return interpolated basises", () => {
        const finalBasis = new Basis(H2HFinal.basis);
        const initialBasis = new Basis(H2HInitial.basis);
        const images = tools.basis.interpolate(initialBasis, finalBasis, 1);
        assertDeepAlmostEqual(H2HImage.basis, images[0].toJSON());
    });
});
