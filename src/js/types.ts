import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";

import { Basis } from "./basis/basis";
import { Lattice } from "./lattice/lattice";
import { MaterialMixinProps } from "./materialMixin";

export type MaterialInMemoryEntity = MaterialMixinProps & InMemoryEntity;
export { Basis };
export { Lattice };
