
//package burai.app.project.viewer.modeler.slabmodel;

//import burai.atoms.model.Cell;

export class SlabModel {

    static get DEFAULT_OFFSET() { return 0.0; }

    static get DEFAULT_THICKNESS() { return 1.0; }

    static get DEFAULT_VACUUM() { return 10.0; } // angstrom

    static get DEFAULT_SCALE() { return 1; }

    static defaultOffset() {
        return SlabModel.DEFAULT_OFFSET;
    }

    static defaultThickness() {
        return SlabModel.DEFAULT_THICKNESS;
    }

    static defaultVacuum() {
        return SlabModel.DEFAULT_VACUUM;
    }

    static defaultScale() {
        return SlabModel.DEFAULT_SCALE;
    }

    #offset = null;
    #thickness = null;
    #vacuum = null;
    #scaleA = null;
    #scaleB = null;

    #lastOffset = null;
    #lastThickness = null;
    #lastVacuum = null;
    #lastScaleA = null;
    #lastScaleB = null;

    constructor() {
        this.offset = SlabModel.DEFAULT_OFFSET;
        this.thickness = SlabModel.DEFAULT_THICKNESS;
        this.vacuum = SlabModel.DEFAULT_VACUUM;
        this.scaleA = SlabModel.DEFAULT_SCALE;
        this.scaleB = SlabModel.DEFAULT_SCALE;

        this.lastOffset = this.offset;
        this.lastThickness = this.thickness;
        this.lastVacuum = this.vacuum;
        this.lastScaleA = this.scaleA;
        this.lastScaleB = this.scaleB;
    }

    setOffset(offset) {
        this.offset = offset;
    }

    getOffset() {
        return this.offset;
    }

    setThickness(thickness) {
        this.thickness = thickness;
    }

    getThickness() {
        return this.thickness;
    }

    setVacuum(vacuum) {
        this.vacuum = vacuum;
    }

    getVacuum() {
        return this.vacuum;
    }

    setScaleA(scaleA) {
        this.scaleA = scaleA;
    }

    getScaleA() {
        return this.scaleA;
    }

    setScaleB(scaleB) {
        this.scaleB = scaleB;
    }

    getScaleB() {
        return this.scaleB;
    }

    //public abstract SlabModel[] getSlabModels();

    //protected abstract boolean updateCell(Cell cell);

    putOnCell(cell) {
        const status = this.updateCell(cell);

        if (status) {
            this.lastOffset = this.offset;
            this.lastThickness = this.thickness;
            this.lastVacuum = this.vacuum;
            this.lastScaleA = this.scaleA;
            this.lastScaleB = this.scaleB;
        }

        return status;
    }

    putOnLastCell(cell) {
        this.offset = this.lastOffset;
        this.thickness = this.lastThickness;
        this.vacuum = this.lastVacuum;
        this.scaleA = this.lastScaleA;
        this.scaleB = this.lastScaleB;

        return this.updateCell(cell);
    }
}
