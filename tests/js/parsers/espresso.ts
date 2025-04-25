import { MaterialSchema } from "@mat3ra/esse/dist/js/types";
import { expect } from "chai";

import { Material } from "../../../src/js/material";
import parsers from "../../../src/js/parsers/parsers";
import { Silicon, SiPWSCFInput } from "../fixtures";

describe("Parsers:Espresso", () => {
    it("should return textual representation of a material according to QE pw.x input format", () => {
        const material = new Material(Silicon);
        expect(parsers.espresso.toEspressoFormat(material as MaterialSchema)).to.be.equal(
            SiPWSCFInput,
        );
    });
});
