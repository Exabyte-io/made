// TODO: refactor and remove the need for UnitCell
export class UnitCell {
    constructor(ax, ay, az, bx, by, bz, cx, cy, cz, units) {
        this.ax = ax;
        this.ay = ay;
        this.az = az;
        this.bx = bx;
        this.by = by;
        this.bz = bz;
        this.cx = cx;
        this.cy = cy;
        this.cz = cz;
        this.units = units;
    }

    vectorA() {
        return [this.ax, this.ay, this.az];
    }

    vectorB() {
        return [this.bx, this.by, this.bz];
    }

    vectorC() {
        return [this.cx, this.cy, this.cz];
    }

    axes() {
        return [this.vectorA(), this.vectorB(), this.vectorC()];
    }
}
