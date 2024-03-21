import { Vector } from "./types";

export type UnitCellProps = [
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    number,
    string,
];

// TODO: refactor and remove the need for UnitCell
export class UnitCell {
    ax: number;

    ay: number;

    az: number;

    bx: number;

    by: number;

    bz: number;

    cx: number;

    cy: number;

    cz: number;

    units: string;

    constructor([ax, ay, az, bx, by, bz, cx, cy, cz, units]: UnitCellProps) {
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

    vectorA(): Vector {
        return [this.ax, this.ay, this.az];
    }

    vectorB(): Vector {
        return [this.bx, this.by, this.bz];
    }

    vectorC(): Vector {
        return [this.cx, this.cy, this.cz];
    }

    axes(): [Vector, Vector, Vector] {
        return [this.vectorA(), this.vectorB(), this.vectorC()];
    }
}
