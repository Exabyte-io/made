import {expect} from "chai";

import {Si, SiPWSCFInput} from "../enums";
import {Material} from "../../src/material";
import parsers from "../../src/parsers/parsers";

describe('Parsers:Espresso', function () {

    it('should return textual representation of a material according to QE pw.x input format', function () {
        const material = new Material(Si);
        expect(parsers.espresso.toEspressoFormat(material)).to.be.equal(SiPWSCFInput)
    });

});
