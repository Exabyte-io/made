import chai from "chai";
import chaiAlmost from "chai-almost";

import { TOLERANCE } from "./enums";

chai.use(chaiAlmost(TOLERANCE));
