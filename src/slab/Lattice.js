
//package burai.com.math;

//import burai.com.consts.Constants;

import {Matrix3D} from './Matrix3D';
import {Constants} from './Constants';

export class Lattice {

    static get ROOT2() { return Math.sqrt(2.0); }
    static get ROOT3() { return Math.sqrt(3.0); }

    static get CELL_THRESHOLD() { return 1.0e-6; }

    //private static final int[] IBRAV_LIST = { 1, 2, 3, -3, 4, 5, -5, 6, 7, 8, 9, -9, 10, 11, 12, -12, 13, 14 };
    static get IBRAV_LIST() { return [ 1, 2, 3, 4, 5, -5, 6, 7, 8, 9, -9, 10, 11, 12, -12, 13, 14 ]; }

    constructor() {
        // NOP
    }

    static checkCell(cell) {
        if (cell === null || cell.length < 3) {
            return false;
        }

        if (cell[0] === null || cell[0].length < 3) {
            return false;
        }

        if (cell[1] === null || cell[1].length < 3) {
            return false;
        }

        if (cell[2] === null || cell[2].length < 3) {
            return false;
        }

        return true;
    }

    static getA(cell) {
        if (!Lattice.checkCell(cell)) {
            return -1.0;
        }

        return Matrix3D.norm(cell[0]);
    }

    static getB(cell) {
        if (!Lattice.checkCell(cell)) {
            return -1.0;
        }

        return Matrix3D.norm(cell[1]);
    }

    static getC(cell) {
        if (!Lattice.checkCell(cell)) {
            return -1.0;
        }

        return Matrix3D.norm(cell[2]);
    }


    static getCosAlpha(cell) {
        if (!Lattice.checkCell(cell)) {
            return 1.0;
        }

        const b = Matrix3D.norm(cell[1]);
        if (b <= 0.0) {
            return 1.0;
        }

        const c = Matrix3D.norm(cell[2]);
        if (c <= 0.0) {
            return 1.0;
        }

        return Matrix3D.mult1D1D(cell[1], cell[2]) / b / c;
    }
/*
    public static double getAlpha(double[][] cell) {
        double cosbc = getCosAlpha(cell);
        return Math.acos(Math.max(-1.0, Math.min(cosbc, 1.0))) * 180.0 / Math.PI;
    }
*/

    static getCosBeta(cell) {
        if (!Lattice.checkCell(cell)) {
            return 1.0;
        }

        const a = Matrix3D.norm(cell[0]);
        if (a <= 0.0) {
            return 1.0;
        }

        const c = Matrix3D.norm(cell[2]);
        if (c <= 0.0) {
            return 1.0;
        }

        return Matrix3D.mult1D1D(cell[0], cell[2]) / a / c;
    }
/*
    public static double getBeta(double[][] cell) {
        double cosac = getCosBeta(cell);
        return Math.acos(Math.max(-1.0, Math.min(cosac, 1.0))) * 180.0 / Math.PI;
    }
*/
    static getCosGamma(cell) {
        if (!Lattice.checkCell(cell)) {
            return 1.0;
        }

        const a = Matrix3D.norm(cell[0]);
        if (a <= 0.0) {
            return 1.0;
        }

        const b = Matrix3D.norm(cell[1]);
        if (b <= 0.0) {
            return 1.0;
        }

        return Matrix3D.mult1D1D(cell[0], cell[1]) / a / b;
    }
/*
    public static double getGamma(double[][] cell) {
        double cosab = getCosGamma(cell);
        return Math.acos(Math.max(-1.0, Math.min(cosab, 1.0))) * 180.0 / Math.PI;
    }

    public static double getXMax(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double x = 0.0;

        if (cell[0][0] > 0.0) {
            x += cell[0][0];
        }

        if (cell[1][0] > 0.0) {
            x += cell[1][0];
        }

        if (cell[2][0] > 0.0) {
            x += cell[2][0];
        }

        return x;
    }

    public static double getXMin(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double x = 0.0;

        if (cell[0][0] < 0.0) {
            x += cell[0][0];
        }

        if (cell[1][0] < 0.0) {
            x += cell[1][0];
        }

        if (cell[2][0] < 0.0) {
            x += cell[2][0];
        }

        return x;
    }

    public static double getYMax(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double y = 0.0;

        if (cell[0][1] > 0.0) {
            y += cell[0][1];
        }

        if (cell[1][1] > 0.0) {
            y += cell[1][1];
        }

        if (cell[2][1] > 0.0) {
            y += cell[2][1];
        }

        return y;
    }

    public static double getYMin(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double y = 0.0;

        if (cell[0][1] < 0.0) {
            y += cell[0][1];
        }

        if (cell[1][1] < 0.0) {
            y += cell[1][1];
        }

        if (cell[2][1] < 0.0) {
            y += cell[2][1];
        }

        return y;
    }

    public static double getZMax(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double z = 0.0;

        if (cell[0][2] > 0.0) {
            z += cell[0][2];
        }

        if (cell[1][2] > 0.0) {
            z += cell[1][2];
        }

        if (cell[2][2] > 0.0) {
            z += cell[2][2];
        }

        return z;
    }

    public static double getZMin(double[][] cell) {
        if (!checkCell(cell)) {
            return 0.0;
        }

        double z = 0.0;

        if (cell[0][2] < 0.0) {
            z += cell[0][2];
        }

        if (cell[1][2] < 0.0) {
            z += cell[1][2];
        }

        if (cell[2][2] < 0.0) {
            z += cell[2][2];
        }

        return z;
    }
*/
    /**
     * celldm of primitive cell
     */
/*
    public static double[] getCellDm(double[][] cell) {
        return getCellDm(0, cell);
    }
*/

    /**
     * celldm of primitive cell (ibrav = 0) or standard cell (ibrav > 0)
     */
    static getCellDm(ibrav, cell) {
        if (!Lattice.checkCell(cell)) {
            return null;
        }

        const celldm = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        const a = Lattice.getA(cell);
        const b = Lattice.getB(cell);
        const c = Lattice.getC(cell);
        const cosAlpha = Lattice.getCosAlpha(cell);
        const cosBeta = Lattice.getCosBeta(cell);
        const cosGamma = Lattice.getCosGamma(cell);
        celldm[0] = a / Constants.BOHR_RADIUS_ANGS;
        celldm[1] = b / a;
        celldm[2] = c / a;
        celldm[3] = cosAlpha;
        celldm[4] = cosBeta;
        celldm[5] = cosGamma;

        if (Lattice.isCorrectBravais(ibrav)) {
            return Lattice.convertCellDm(ibrav, celldm);
        } else {
            return celldm;
        }
    }

    /**
     * a, b, c, alpha, beta, gamma of primitive cell
     */
/*
    public static double[] getLatticeConstants(double[][] cell, boolean asCos) {
        return getLatticeConstants(0, cell, asCos);
    }
*/

    /**
     * a, b, c, alpha, beta, gamma of primitive cell (ibrav = 0) or standard cell (ibrav > 0)
     */
    static getLatticeConstants(ibrav, cell, asCos) {
        if (!Lattice.checkCell(cell)) {
            return null;
        }

        const celldm = Lattice.getCellDm(ibrav, cell);
        if (celldm === null || celldm.length < 6) {
            return null;
        }

        const a = celldm[0] * Constants.BOHR_RADIUS_ANGS;
        const b = a * celldm[1];
        const c = a * celldm[2];

        let cosAlpha = 0.0;
        let cosBeta = 0.0;
        let cosGamma = 0.0;
        if (ibrav == 14) {
            cosAlpha = celldm[3];
            cosBeta = celldm[4];
            cosGamma = celldm[5];
        } else if (ibrav == -12 || ibrav == -13) {
            cosAlpha = 0.0;
            cosBeta = celldm[4];
            cosGamma = 0.0;
        } else if (Lattice.isCorrectBravais(ibrav)) {
            cosAlpha = 0.0;
            cosBeta = 0.0;
            cosGamma = celldm[3];
        } else {
            cosAlpha = celldm[3];
            cosBeta = celldm[4];
            cosGamma = celldm[5];
        }

        if (asCos) {
            return [ a, b, c, cosAlpha, cosBeta, cosGamma ];
        }

        const alpha = Math.acos(Math.max(-1.0, Math.min(cosAlpha, 1.0))) * 180.0 / Math.PI;
        const beta = Math.acos(Math.max(-1.0, Math.min(cosBeta, 1.0))) * 180.0 / Math.PI;
        const gamma = Math.acos(Math.max(-1.0, Math.min(cosGamma, 1.0))) * 180.0 / Math.PI;

        return [ a, b, c, alpha, beta, gamma ];
    }

    /**
     * check available value of ibrav
     */
    static isCorrectBravais(ibrav) {
        for (let ibrav2 of Lattice.IBRAV_LIST) {
            if (ibrav === ibrav2) {
                return true;
            }
        }

        return false;
    }

    /**
     * detect ibrav from lattice vectors
     */
    static getBravais(cell) {
        if (!Lattice.checkCell(cell)) {
            return 0;
        }

        // lattice vector -> primitive celldm
        const celldmPrim = Lattice.getCellDm(0, cell);
        if (celldmPrim === null || celldmPrim.length < 6) {
            return 0;
        }

        for (const ibrav of Lattice.IBRAV_LIST) {
            // primitive celldm -> standard celldm
            const celldmStd = Lattice.convertCellDm(ibrav, celldmPrim);
            if (celldmStd === null || celldmStd.length < 6) {
                continue;
            }

            // standard celldm -> lattice vectors
            const cell_ = Lattice.getCell(ibrav, celldmStd);
            if (cell_ !== null && Matrix3D.equals2D2DThreshold(cell, cell_, Lattice.CELL_THRESHOLD)) {
                return ibrav;
            }
        }

        return 0;
    }

    /**
     * detect lattice vectors from primitive celldm
     */
/*
    private static double[][] getCell(double[] celldmPrim) {
        // primitive celldm
        if (celldmPrim == null || celldmPrim.length < 6) {
            return null;
        }

        for (int ibrav : IBRAV_LIST) {
            // primitive celldm -> standard celldm
            double[] celldmStd = convertCellDm(ibrav, celldmPrim);
            if (celldmStd == null || celldmStd.length < 6) {
                continue;
            }

            // standard celldm -> lattice vectors
            double[][] cell = getCell(ibrav, celldmStd);
            if (cell == null) {
                continue;
            }

            // lattice vectors -> primitive celldm
            double[] celldmPrim_ = getCellDm(cell);
            if (celldmPrim_ == null || celldmPrim_.length < 6) {
                continue;
            }

            boolean sameCell = true;
            for (int i = 0; i < 6; i++) {
                if (Math.abs(celldmPrim[i] - celldmPrim_[i]) > CELL_THRESHOLD) {
                    sameCell = false;
                    break;
                }
            }

            if (sameCell) {
                return cell;
            }
        }

        return null;
    }
*/

    /**
     * detect lattice vectors from primitive a, b, c, alpha, beta, gamma
     */
/*
    public static double[][] getCell(double a, double b, double c, double alpha, double beta, double gamma) {
        if (a <= 0.0) {
            return null;
        }
        if (b <= 0.0) {
            return null;
        }
        if (c <= 0.0) {
            return null;
        }
        if (alpha <= 0.0 || 180.0 <= alpha) {
            return null;
        }
        if (beta <= 0.0 || 180.0 <= beta) {
            return null;
        }
        if (gamma <= 0.0 || 180.0 <= gamma) {
            return null;
        }

        // primitive celldm
        double[] celldm = new double[6];
        double cosAlpha = Math.cos(alpha * Math.PI / 180.0);
        double cosBeta = Math.cos(beta * Math.PI / 180.0);
        double cosGamma = Math.cos(gamma * Math.PI / 180.0);
        celldm[0] = a / Constants.BOHR_RADIUS_ANGS;
        celldm[1] = b / a;
        celldm[2] = c / a;
        celldm[3] = cosAlpha;
        celldm[4] = cosBeta;
        celldm[5] = cosGamma;

        // primitive celldm -> lattice vectors
        return getCell(celldm);
    }
*/

    /**
     * convert celldm: primitive cell -> standard cell
     */
    static convertCellDm(ibrav, celldmPrim) {
        if (celldmPrim === null || celldmPrim.length < 6) {
            return null;
        }

        let celldmStd = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        //for (let i = 0; i < 6; i++) {
        //    celldmStd[i] = 0.0;
        //}

        switch (ibrav) {
        case 1:
            celldmStd[0] = celldmPrim[0];
            break;

        case 2:
            celldmStd[0] = (2.0 / Lattice.ROOT2) * celldmPrim[0];
            break;

        case 3:
            celldmStd[0] = (2.0 / Lattice.ROOT3) * celldmPrim[0];
            break;

        case -3:
            celldmStd[0] = (2.0 / Lattice.ROOT3) * celldmPrim[0];
            break;

        case 4:
            celldmStd[0] = celldmPrim[0];
            celldmStd[2] = celldmPrim[2];
            break;

        case 5:
            celldmStd[0] = celldmPrim[0];
            celldmStd[3] = celldmPrim[5];
            break;

        case -5:
            celldmStd[0] = celldmPrim[0];
            celldmStd[3] = celldmPrim[5];
            break;

        case 6:
            celldmStd[0] = celldmPrim[0];
            celldmStd[2] = celldmPrim[2];
            break;

        case 7:
            if (celldmPrim[0] <= 0.0) {
                celldmStd = null;
                break;
            }
            if (celldmPrim[5] <= 0.0 || 1.0 <= celldmPrim[5]) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = Lattice.ROOT2 * celldmPrim[0] * Math.sqrt(1.0 - celldmPrim[5]);
            celldmStd[2] = 2.0 * celldmPrim[0] * Math.sqrt(celldmPrim[5]) / celldmStd[0];
            break;

        case 8:
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = celldmPrim[1];
            celldmStd[2] = celldmPrim[2];
            break;

        case 9:
            if (celldmPrim[0] <= 0.0) {
                celldmStd = null;
                break;
            }
            if (Math.abs(celldmPrim[5]) >= 1.0) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = Lattice.ROOT2 * celldmPrim[0] * Math.sqrt(1 - celldmPrim[5]);
            celldmStd[1] = Math.sqrt(4.0 * celldmPrim[0] * celldmPrim[0] - celldmStd[0] * celldmStd[0]) / celldmStd[0];
            celldmStd[2] = celldmPrim[2] * celldmPrim[0] / celldmStd[0];
            break;

        case -9:
            if (celldmPrim[0] <= 0.0) {
                celldmStd = null;
                break;
            }
            if (Math.abs(celldmPrim[5]) >= 1.0) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = Lattice.ROOT2 * celldmPrim[0] * Math.sqrt(1 + celldmPrim[5]);
            celldmStd[1] = Math.sqrt(4.0 * celldmPrim[0] * celldmPrim[0] - celldmStd[0] * celldmStd[0]) / celldmStd[0];
            celldmStd[2] = celldmPrim[2] * celldmPrim[0] / celldmStd[0];
            break;

        case 91:
            if (Math.abs(celldmPrim[3]) >= 1.0) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = Lattice.ROOT2 * celldmPrim[1] * Math.sqrt(1 + celldmPrim[3]);
            celldmStd[2] = Math.sqrt(4.0 * celldmPrim[1] * celldmPrim[1] - celldmStd[1] * celldmStd[1]);
            break;

        case 10:
            const ap = celldmPrim[0];
            const bp = celldmPrim[1] * ap;
            const cp = celldmPrim[2] * ap;
            const asSqr = ap * ap + bp * bp - cp * cp;
            const bsSqr = bp * bp + cp * cp - ap * ap;
            const csSqr = ap * ap + cp * cp - bp * bp;
            if (asSqr <= 0.0 || bsSqr <= 0.0 || csSqr <= 0.0) {
                celldmStd = null;
                break;
            }

            const as = Lattice.ROOT2 * Math.sqrt(asSqr);
            const bs = Lattice.ROOT2 * Math.sqrt(bsSqr);
            const cs = Lattice.ROOT2 * Math.sqrt(csSqr);
            celldmStd[0] = as;
            celldmStd[1] = bs / as;
            celldmStd[2] = cs / as;
            break;

        case 11:
            const rr = 4.0 * celldmPrim[0] * celldmPrim[0];
            const x1 = rr * celldmPrim[3];
            const x2 = rr * celldmPrim[4];
            const x3 = rr * celldmPrim[5];
            const aa = 0.5 * (x1 - x2);
            const bb = 0.5 * (x2 - x3);
            const cc = 0.5 * (x3 + x1);
            if (aa <= 0.0 || bb <= 0.0 || cc <= 0.0) {
                celldmStd = null;
                break;
            }

            const a = Math.sqrt(aa);
            const b = Math.sqrt(bb);
            const c = Math.sqrt(cc);
            celldmStd[0] = a;
            celldmStd[1] = b / a;
            celldmStd[2] = c / a;
            break;

        case 12:
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = celldmPrim[1];
            celldmStd[2] = celldmPrim[2];
            celldmStd[3] = celldmPrim[5];
            break;

        case -12:
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = celldmPrim[1];
            celldmStd[2] = celldmPrim[2];
            celldmStd[4] = celldmPrim[4];
            break;

        case 13:
            if (celldmPrim[0] <= 0.0) {
                celldmStd = null;
                break;
            }
            if (Math.abs(celldmPrim[4]) >= 1.0) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = Lattice.ROOT2 * celldmPrim[0] * Math.sqrt(1 + celldmPrim[4]);
            celldmStd[1] = celldmPrim[1] * celldmPrim[0] / celldmStd[0];
            celldmStd[2] = Math.sqrt(4.0 * celldmPrim[0] * celldmPrim[0] - celldmStd[0] * celldmStd[0]) / celldmStd[0];
            celldmStd[3] = celldmPrim[5] * Math.sqrt(1.0 + celldmStd[2] * celldmStd[2]);
            break;

        case -13:
            if (celldmPrim[0] <= 0.0) {
                celldmStd = null;
                break;
            }
            if (Math.abs(celldmPrim[5]) >= 1.0) {
                celldmStd = null;
                break;
            }
            celldmStd[0] = Lattice.ROOT2 * celldmPrim[0] * Math.sqrt(1 + celldmPrim[5]);
            celldmStd[1] = Math.sqrt(4.0 * celldmPrim[0] * celldmPrim[0] - celldmStd[0] * celldmStd[0]) / celldmStd[0];
            celldmStd[2] = celldmPrim[2] * celldmPrim[0] / celldmStd[0];
            celldmStd[4] = celldmPrim[4] * Math.sqrt(1.0 + celldmStd[1] * celldmStd[1]);
            break;

        case 14:
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = celldmPrim[1];
            celldmStd[2] = celldmPrim[2];
            celldmStd[3] = celldmPrim[3];
            celldmStd[4] = celldmPrim[4];
            celldmStd[5] = celldmPrim[5];
            break;

        default:
            celldmStd[0] = celldmPrim[0];
            celldmStd[1] = celldmPrim[1];
            celldmStd[2] = celldmPrim[2];
            celldmStd[3] = celldmPrim[3];
            celldmStd[4] = celldmPrim[4];
            celldmStd[5] = celldmPrim[5];
            break;
        }

        return celldmStd;
    }

    /**
     * create lattice vectors from standard celldm
     */
    static getCell(ibrav, celldm) {
        if (celldm === null || celldm.length < 6) {
            return null;
        }

        if (celldm[0] === 0.0) {
            return null;
        }

        let lattice = Matrix3D.zero();

        let term1;
        let term2;

        switch (ibrav) {
        case 1:
            lattice[0][0] = celldm[0];
            lattice[1][1] = celldm[0];
            lattice[2][2] = celldm[0];
            break;

        case 2:
            term1 = 0.5 * celldm[0];
            lattice[0][0] = -term1;
            lattice[0][2] = term1;
            lattice[1][1] = term1;
            lattice[1][2] = term1;
            lattice[2][0] = -term1;
            lattice[2][1] = term1;
            break;

        case 3:
            term1 = 0.5 * celldm[0];
            for (let i = 0; i < 3; i++) {
                lattice[0][i] = term1;
                lattice[1][i] = term1;
                lattice[2][i] = term1;
            }
            lattice[1][0] *= -1.0;
            lattice[2][0] *= -1.0;
            lattice[2][1] *= -1.0;
            break;

        case -3:
            term1 = 0.5 * celldm[0];
            for (let i = 0; i < 3; i++) {
                lattice[0][i] = term1;
                lattice[1][i] = term1;
                lattice[2][i] = term1;
            }
            lattice[0][0] *= -1.0;
            lattice[1][1] *= -1.0;
            lattice[2][2] *= -1.0;
            break;

        case 4:
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][0] = -celldm[0] / 2.0;
            lattice[1][1] = celldm[0] * Math.sqrt(3.0) / 2.0;
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case 5:
            if (celldm[3] <= -0.5 || celldm[3] >= 1.0) {
                lattice = null;
                break;
            }
            term1 = Math.sqrt(1.0 + 2.0 * celldm[3]);
            term2 = Math.sqrt(1.0 - celldm[3]);
            lattice[1][1] = Lattice.ROOT2 * celldm[0] * term2 / Lattice.ROOT3;
            lattice[1][2] = celldm[0] * term1 / Lattice.ROOT3;
            lattice[0][0] = celldm[0] * term2 / Lattice.ROOT2;
            lattice[0][1] = -lattice[0][0] / Lattice.ROOT3;
            lattice[0][2] = lattice[1][2];
            lattice[2][0] = -lattice[0][0];
            lattice[2][1] = lattice[0][1];
            lattice[2][2] = lattice[1][2];
            break;

        case -5:
            if (celldm[3] <= -0.5 || celldm[3] >= 1.0) {
                lattice = null;
                break;
            }
            term1 = Math.sqrt(1.0 + 2.0 * celldm[3]);
            term2 = Math.sqrt(1.0 - celldm[3]);
            lattice[0][0] = celldm[0] * (term1 - 2.0 * term2) / 3.0;
            lattice[0][1] = celldm[0] * (term1 + term2) / 3.0;
            lattice[0][2] = lattice[0][1];
            lattice[1][0] = lattice[0][2];
            lattice[1][1] = lattice[0][0];
            lattice[1][2] = lattice[0][1];
            lattice[2][0] = lattice[0][1];
            lattice[2][1] = lattice[0][2];
            lattice[2][2] = lattice[0][0];
            break;

        case 6:
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][1] = celldm[0];
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case 7:
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[1][0] = celldm[0] / 2.0;
            lattice[1][1] = lattice[1][0];
            lattice[1][2] = celldm[2] * celldm[0] / 2.0;
            lattice[0][0] = lattice[1][0];
            lattice[0][1] = -lattice[1][0];
            lattice[0][2] = lattice[1][2];
            lattice[2][0] = -lattice[1][0];
            lattice[2][1] = -lattice[1][0];
            lattice[2][2] = lattice[1][2];
            break;

        case 8:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][1] = celldm[0] * celldm[1];
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case 9:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = 0.5 * celldm[0];
            lattice[0][1] = lattice[0][0] * celldm[1];
            lattice[1][0] = -lattice[0][0];
            lattice[1][1] = lattice[0][1];
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case -9:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = 0.5 * celldm[0];
            lattice[0][1] = -lattice[0][0] * celldm[1];
            lattice[1][0] = lattice[0][0];
            lattice[1][1] = -lattice[0][1];
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case 91:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][1] = celldm[0] * celldm[1] * 0.5;
            lattice[1][2] = -celldm[0] * celldm[2] * 0.5;
            lattice[2][1] = lattice[1][1];
            lattice[2][2] = -lattice[1][2];
            break;

        case 10:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[1][0] = 0.5 * celldm[0];
            lattice[1][1] = lattice[1][0] * celldm[1];
            lattice[0][0] = lattice[1][0];
            lattice[0][2] = lattice[1][0] * celldm[2];
            lattice[2][1] = lattice[1][0] * celldm[1];
            lattice[2][2] = lattice[0][2];
            break;

        case 11:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = 0.5 * celldm[0];
            lattice[0][1] = lattice[0][0] * celldm[1];
            lattice[0][2] = lattice[0][0] * celldm[2];
            lattice[1][0] = -lattice[0][0];
            lattice[1][1] = lattice[0][1];
            lattice[1][2] = lattice[0][2];
            lattice[2][0] = -lattice[0][0];
            lattice[2][1] = -lattice[0][1];
            lattice[2][2] = lattice[0][2];
            break;

        case 12:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[3]) >= 1.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][0] = celldm[0] * celldm[1] * celldm[3];
            lattice[1][1] = celldm[0] * celldm[1] * Math.sqrt(1.0 - celldm[3] * celldm[3]);
            lattice[2][2] = celldm[0] * celldm[2];
            break;

        case -12:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[4]) >= 1.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = celldm[0];
            lattice[1][1] = celldm[0] * celldm[1];
            lattice[2][0] = celldm[0] * celldm[2] * celldm[4];
            lattice[2][2] = celldm[0] * celldm[2] * Math.sqrt(1.0 - celldm[4] * celldm[4]);
            break;

        case 13:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[3]) >= 1.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = 0.5 * celldm[0];
            lattice[0][2] = -lattice[0][0] * celldm[2];
            lattice[1][0] = celldm[0] * celldm[1] * celldm[3];
            lattice[1][1] = celldm[0] * celldm[1] * Math.sqrt(1.0 - celldm[3] * celldm[3]);
            lattice[2][0] = lattice[0][0];
            lattice[2][2] = -lattice[0][2];
            break;

        case -13:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[4]) >= 1.0) {
                lattice = null;
                break;
            }
            lattice[0][0] = 0.5 * celldm[0];
            lattice[0][1] = -lattice[0][0] * celldm[1];
            lattice[1][0] = lattice[0][0];
            lattice[1][1] = -lattice[0][1];
            lattice[2][0] = celldm[0] * celldm[2] * celldm[4];
            lattice[2][2] = celldm[0] * celldm[2] * Math.sqrt(1.0 - celldm[4] * celldm[4]);
            break;

        case 14:
            if (celldm[1] <= 0.0) {
                lattice = null;
                break;
            }
            if (celldm[2] <= 0.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[3]) >= 1.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[4]) >= 1.0) {
                lattice = null;
                break;
            }
            if (Math.abs(celldm[5]) >= 1.0) {
                lattice = null;
                break;
            }
            term1 = Math.sqrt(1.0 - celldm[5] * celldm[5]);
            if (term1 === 0.0) {
                lattice = null;
                break;
            }
            let term2 = 1.0 + 2.0 * celldm[3] * celldm[4] * celldm[5];
            term2 += -celldm[3] * celldm[3] - celldm[4] * celldm[4] - celldm[5] * celldm[5];
            if (term2 < 0.0) {
                lattice = null;
                break;
            }
            term2 = Math.sqrt(term2 / (1.0 - celldm[5] * celldm[5]));
            lattice[0][0] = celldm[0];
            lattice[1][0] = celldm[0] * celldm[1] * celldm[5];
            lattice[1][1] = celldm[0] * celldm[1] * term1;
            lattice[2][0] = celldm[0] * celldm[2] * celldm[4];
            lattice[2][1] = celldm[0] * celldm[2] * (celldm[3] - celldm[4] * celldm[5]) / term1;
            lattice[2][2] = celldm[0] * celldm[2] * term2;
            break;

        default:
            lattice = null;
            break;
        }

        if (lattice != null) {
            lattice = Matrix3D.multAlpha2D(Constants.BOHR_RADIUS_ANGS, lattice);
        }

        return lattice;
    }
}
