import { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import { MaterialSchema } from "@mat3ra/esse/dist/js/types";

import { Basis } from "../basis/basis";
import { ConstrainedBasis } from "../basis/constrained_basis";
import { Cell } from "../cell/cell";

export type MaterialJSON = MaterialSchema & AnyObject;

/**
 * Interface defining the minimum required properties for a Material
 * when used by other components like the supercell tools.
 * This helps break circular dependencies.
 */
export interface MaterialInterface {
    name: string;
    Basis: Basis | ConstrainedBasis;
    Lattice: {
        vectors: Cell;
    };
}
