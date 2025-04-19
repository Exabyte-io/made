import chai from "chai";
import chaiAlmost from "chai-almost";

const TOLERANCE = 1e-3;

chai.use(chaiAlmost(TOLERANCE));
