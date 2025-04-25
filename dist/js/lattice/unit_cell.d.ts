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
    string
];
export declare class UnitCell {
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
    constructor([ax, ay, az, bx, by, bz, cx, cy, cz, units]: UnitCellProps);
    vectorA(): Vector3DSchema;
    vectorB(): Vector3DSchema;
    vectorC(): Vector3DSchema;
    axes(): [Vector3DSchema, Vector3DSchema, Vector3DSchema];
}
