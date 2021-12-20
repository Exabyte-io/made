
//package burai.app.project.viewer.modeler.slabmodel;

//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.List;

export class SlabGenom {

    static get COORD_THR() { return 0.10; } // angstrom
    static get LAYER_THR() { return 0.10; } // angstrom

    #layers = null;

    constructor(names, coords) {
        if (names === null || names.length < 1) {
            throw new IllegalArgumentException("names is empty.");
        }

        if (coords === null || coords.length < 1) {
            throw new IllegalArgumentException("coords is empty.");
        }

        if (names.length !== coords.length) {
            throw new IllegalArgumentException("names.length != coords.length.");
        }

        this.setupLayers(names, coords);
    }

    setupLayers(names, coords) {
        this.layers = [];

        let istart = 0;
        let iend = 0;
        let coord1 = 0.0;
        let coord2 = 0.0;

        while (true) {
            istart = iend;
            iend = this.nextLayer(istart, coords);
            if (istart >= iend) {
                break;
            }

            coord1 = coord2;
            coord2 = 0.0;
            for (let i = istart; i < iend; i++) {
                coord2 += coords[i];
            }
            coord2 /= iend - istart;

            let distance = 0.0;
            if (this.layers.length !== 0) {
                distance = coord1 - coord2;
            }

            const layer = this.getLayer(istart, iend, names);
            if (layer !== null) {
                layer.distance = distance;
                this.layers.push(layer);
            }
        }
    }

    nextLayer(istart, coords) {
        if (istart >= coords.length) {
            return coords.length;
        }

        let iend = coords.length;
        const coord0 = coords[istart];

        for (let i = istart + 1; i < coords.length; i++) {
            const coord = coords[i];
            if (Math.abs(coord - coord0) > SlabGenom.COORD_THR) {
                iend = i;
                break;
            }
        }

        return iend;
    }

    getLayer(istart, iend, names) {
        if ((iend - istart) === 1) {
            const layer = new Layer();
            layer.code = names[istart];
            return layer;
        }

        const names2 = [];
        for(let i = 0; i < iend - istart; i++) names2.push('');
        for (let i = istart; i < iend; i++) {
            names2[i - istart] = names[i];
        }

        names2.sort();

        let mult = 0;
        let name = null;
        let code = '';

        for (let i = 0; i <= names2.length; i++) {
            if (i < names2.length && name != null && name === names2[i]) {
                mult++;
                continue;
            }

            if (name !== null && !(name.length === 0)) {
                if (code.length > 0) {
                    code += ' ';
                }

                code += name;
                if (mult > 1) {
                    code += '*';
                    code += mult;
                }
            }

            if (i < names2.length) {
                mult = 1;
                name = names2[i];
            }
        }

        const layer = new Layer();
        layer.code = code;
        return layer;
    }

/*
    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();

        if (this.layers != null) {
            for (Layer layer : this.layers) {
                if (layer != null) {
                    str.append('{');
                    str.append(layer.code);
                    str.append('|');
                    str.append(layer.distance);
                    str.append('}');
                }
            }
        }

        return str.length() > 0 ? str.toString() : "{}";
    }

    @Override
    public int hashCode() {
        return this.layers == null ? 0 : this.layers.hashCode();
    }

*/

//    @Override
    equals(obj) {
        if (this === obj) {
            return true;
        }

        if (obj === null) {
            return false;
        }

        //if (this.getClass() != obj.getClass()) {
        //    return false;
        //}

        const other = obj;
        if (this.layers === null) {
            return this.layers === other.layers;
        } else {
            if (this.layers.length !== other.layers.length) {
                return false;
            }

            for (const i in this.layers) {
                if (!(this.layers[i].equals(other.layers[i]))) {
                    return false;
                }
            }

            return true;
//            return this.layers.equals(other.layers);
        }
    }
}

class Layer {
    code = null;
    distance = null;

    constructor() {
        // NOP
    }

/*
    @Override
    public int hashCode() {
        return this.code == null ? 0 : this.code.hashCode();
    }
*/

//    @Override
    equals(obj) {
        if (this === obj) {
            return true;
        }

        if (obj === null) {
            return false;
        }

        //if (this.getClass() != obj.getClass()) {
        //    return false;
        //}

        const other = obj;

        if (this.code === null) {
            if (this.code !== other.code) {
                return false;
            }
        } else {
            if (this.code !== other.code) {
                return false;
            }
        }

        return Math.abs(this.distance - other.distance) <= SlabGenom.LAYER_THR;
    }
}

