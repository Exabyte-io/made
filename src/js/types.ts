import { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import { MaterialSchema } from "@mat3ra/esse/dist/js/types";

import { Basis } from "./basis/basis";
import { ConstrainedBasis } from "./basis/constrained_basis";
import { Lattice } from "./lattice/lattice";

export type MaterialJSON = MaterialSchema & AnyObject;

export interface MaterialInterface {
    name: string;
    Basis: Basis | ConstrainedBasis;
    Lattice: Lattice;
}

/**
 * @deprecated Import from './types' directory instead
 */

export * from './types/material';
