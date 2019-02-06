import {ArrayWithIds} from "../abstract/array_with_ids";

export class AtomicConstraints {

    static fromArray(array) {return new AtomicConstraints({values: array})}

    constructor(config = {}) {
        this.name = 'atomic_constraints';
        this.values = new ArrayWithIds(config.values || []);
    }

    toJSON() {
        return {
            name: this.name,
            values: this.values.toJSON()
        }
    }

    getByIndex(idx) {return this.values.getArrayElementByIndex(idx) || []}

    getAsStringByIndex(idx, func = (boolean) => boolean ? 1 : 0) {
        return this.getByIndex(idx).map(func).join(' ');
    }

}
