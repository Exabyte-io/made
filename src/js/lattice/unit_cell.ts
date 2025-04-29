import { Vector3DSchema } from "@mat3ra/esse/dist/js/types";

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

    vectorA(): Vector3DSchema {
        return [this.ax, this.ay, this.az];
    }

    vectorB(): Vector3DSchema {
        return [this.bx, this.by, this.bz];
    }

    vectorC(): Vector3DSchema {
        return [this.cx, this.cy, this.cz];
    }

    axes(): [Vector3DSchema, Vector3DSchema, Vector3DSchema] {
        return [this.vectorA(), this.vectorB(), this.vectorC()];
    }
}
