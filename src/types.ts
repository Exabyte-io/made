import { AnyObject } from "@exabyte-io/code.js/dist/entity/in_memory";

import { ConstrainedBasisJSON } from "./basis/constrained_basis";
import { LatticeJSON } from "./lattice/lattice";

export interface MaterialJSON extends AnyObject {
    lattice: LatticeJSON;
    basis: ConstrainedBasisJSON;
    name: string;
    isNonPeriodic: boolean;
}
