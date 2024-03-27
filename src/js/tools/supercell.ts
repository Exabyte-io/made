import { Basis } from "../basis/basis";
import { Coordinate } from "../basis/types";
import { Cell } from "../cell/cell";
import { LatticeBravais } from "../lattice/lattice_bravais";
// eslint-disable-next-line import/no-cycle
import { Material } from "../material";
import math from "../math";
import cellTools from "./cell";

const ADD = math.add;

/**
 * @summary Generates new basis for a supercell. For each site from basis generates shifts that are within supercell.
 */
function generateNewBasisWithinSupercell(
    basis: Basis,
    cell: Cell,
    supercell: Cell,
    supercellMatrix: number[][],
): Basis {
    const oldBasis = basis.clone();
    const newBasis = basis.clone({ isEmpty: true });

    oldBasis.toCrystal();
    newBasis.toCrystal();

    oldBasis.elements.forEach((element) => {
        const coordinate = oldBasis.getCoordinateByIndex(element.id);
        const cartesianCoordinate = cell.convertPointToCartesian(coordinate);
        const shifts = cellTools.latticePointsInSupercell(supercellMatrix);
        shifts.forEach((comb) => {
            // "combination" is effectively a point in fractional coordinates here, hence the below
            const newPoint = ADD(cartesianCoordinate, supercell.convertPointToCartesian(comb));
            if (supercell.isPointInsideCell(newPoint)) {
                newBasis.addAtom({
                    element: element.value,
                    coordinate: supercell.convertPointToFractional(newPoint) as Coordinate,
                });
            }
        });
    });

    return newBasis;
}

/**
 * @summary Generates supercell config for the specified material.
 * @param material
 * @param supercellMatrix {Number[][]}
 */
function generateConfig(material: Material, supercellMatrix: number[][]) {
    const det = math.det(supercellMatrix);
    if (det === 0) {
        throw new Error("Scaling matrix is degenerate.");
    }
    const cell = material.Lattice.Cell;
    const supercell = cell.cloneAndScaleByMatrix(supercellMatrix);
    const newBasis = generateNewBasisWithinSupercell(
        material.Basis,
        cell,
        supercell,
        supercellMatrix,
    );
    const newLattice = LatticeBravais.fromVectors({
        a: supercell.vector1,
        b: supercell.vector2,
        c: supercell.vector3,
    });

    return {
        name: `${material.name} - supercell ${JSON.stringify(supercellMatrix)}`,
        basis: newBasis.toJSON(),
        lattice: newLattice.toJSON(),
    };
}

export default {
    generateConfig,
    generateNewBasisWithinSupercell,
};
