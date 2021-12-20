
//package burai.app.project.viewer.modeler.slabmodel;

//import burai.atoms.model.Cell;

import {SlabModel} from './SlabModel';

export class SlabModelLeaf extends SlabModel {

    #stem = null;

    constructor(stem, offset) {
        super();

        if (stem === null) {
            throw new IllegalArgumentException("stem is null.");
        }

        this.stem = stem;
        this.offset = offset;
    }

//    @Override
//    public SlabModel[] getSlabModels() {
//        return new SlabModel[] { this };
//    }

//    @Override
    updateCell(cell) {
        return this.stem.updateCell(cell, this);
    }

}
