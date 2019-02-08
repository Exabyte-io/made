import s from "underscore.string";
import _ from "underscore";
import xyz from "./xyz";
import {Lattice} from "../lattice/lattice";

function toEspressoFormat(material) {
    const l = new Lattice(material.lattice);
    const vectors = l.vectorArrays;

    const vectorsAsString = _.map(vectors, function (v) {
        return s.sprintf('%14.9f', v[0]) + '\t' + s.sprintf('%14.9f', v[1]) + '\t' + s.sprintf('%14.9f', v[2]);
    }).join('\n');
    return s.sprintf('CELL_PARAMETERS (angstroms)\n%s\n\nATOMIC_POSITIONS (crystal)\n%s',
        vectorsAsString, xyz.fromMaterial(material));
}

export default {
    toEspressoFormat
}
