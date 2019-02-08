import s from "underscore.string";

import math from "../math";
import {Lattice} from "../lattice/lattice";
import {ATOMIC_COORD_UNITS} from "../constants";
import {ConstrainedBasis} from "../basis/constrained_basis";

const _print = (x, printFormat = '%14.9f') => s.sprintf(printFormat, math.precise(x));

const _latticeVectorsToString = (vectors) => vectors.map((v) => v.map(c => _print(c)).join('\t')).join('\n');

const _atomicConstraintsCharFromBool = (bool) => bool ? "T" : "F";

function toPoscar(materialConfig, omitConstraints = false) {
    const lattice = new Lattice(materialConfig.lattice);
    const vectorsAsString = _latticeVectorsToString(lattice.vectorArrays);

    const basis = new ConstrainedBasis({
        ...materialConfig.basis,
        cell: lattice.vectorArrays
    });

    const BasisLines = [];
    let addSelectiveDynamics = false;
    basis._elements.array.forEach((item, idx) => {
        const coord = basis.getCoordinateByIndex(idx).map(x => _print(x));
        const constraintsAsString = omitConstraints ?
            '' : basis.AtomicConstraints.getAsStringByIndex(idx, _atomicConstraintsCharFromBool);
        if (constraintsAsString && !omitConstraints) addSelectiveDynamics = true;
        BasisLines.push([coord.join(' '), constraintsAsString, item].join(' '));
    });
    const basisContent = BasisLines.join('\n');

    const elementsLine = Object.keys(basis.elementCounts).join(' ');
    const countsLine = Object.values(basis.elementCounts).map(x => parseInt(x)).join(' ');

    const coordsType = materialConfig.basis.units === ATOMIC_COORD_UNITS.cartesian ? 'cartesian' : 'direct';

    return [
        materialConfig.name,
        "1.0",
        vectorsAsString,
        elementsLine,
        countsLine,
        // add selective dynamics only if there are some constraints!
        ...addSelectiveDynamics ? ["Selective dynamics"] : [],
        coordsType,
        basisContent
    ].join('\n');
}

export default {
    toPoscar,
}
