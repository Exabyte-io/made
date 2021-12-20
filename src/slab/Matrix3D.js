
//package burai.com.math;

export class Matrix3D {

    static get MIN_DET() { return 1.0e-20; }

    static get EQUALS_THRESHOLD() { return 1.0e-10; }

    constructor() {
        // NOP
    }

    static checkMatrix2D(matrix) {
        if (matrix === null || matrix.length < 3) {
            return false;
        }

        for (let i = 0; i < 3; i++) {
            if (matrix[i] === null || matrix[i].length < 3) {
                return false;
            }
        }

        return true;
    }

    static checkMatrix1D(matrix) {
        if (matrix === null || matrix.length < 3) {
            return false;
        }

        return true;
    }

    static zero1() {
        const matrix = [0.0, 0.0, 0.0];
        //for (let i = 0; i < 3; i++) {
        //    matrix[i] = 0.0;
        //}

        return matrix;
    }

    static zero() {
        const matrix = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        //for (let i = 0; i < 3; i++) {
        //    for (let j = 0; j < 3; j++) {
        //        matrix[i][j] = 0.0;
        //    }
        //}

        return matrix;
    }

    static unitOne() {
        return Matrix3D.unit(1.0);
    }

    static unit(alpha) {
        const matrix = Matrix3D.zero();
        matrix[0][0] = alpha;
        matrix[1][1] = alpha;
        matrix[2][2] = alpha;

        return matrix;
    }

    static plus1D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix1D(matrix1) && Matrix3D.checkMatrix1D(matrix2))) {
            return null;
        }

        const matrix3 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix3[i] = matrix1[i] + matrix2[i];
        }

        return matrix3;
    }

    static plus2D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix2D(matrix1) && Matrix3D.checkMatrix2D(matrix2))) {
            return null;
        }

        const matrix3 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix3[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }

        return matrix3;
    }

    static minus1D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix1D(matrix1) && Matrix3D.checkMatrix1D(matrix2))) {
            return null;
        }

        const matrix3 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix3[i] = matrix1[i] - matrix2[i];
        }

        return matrix3;
    }

    static minus2D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix2D(matrix1) && Matrix3D.checkMatrix2D(matrix2))) {
            return null;
        }

        const matrix3 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix3[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }

        return matrix3;
    }

    static trans(matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return null;
        }

        const matrix2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix2[i][j] = matrix[j][i];
            }
        }

        return matrix2;
    }

    static mult1D1D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix1D(matrix1) && Matrix3D.checkMatrix1D(matrix2))) {
            return 0.0;
        }

        let dot = 0.0;
        for (let i = 0; i < 3; i++) {
            dot += matrix1[i] * matrix2[i];
        }

        return dot;
    }

    static mult2D1D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix2D(matrix1) && Matrix3D.checkMatrix1D(matrix2))) {
            return null;
        }

        const matrix3 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix3[i] = 0.0;
            for (let j = 0; j < 3; j++) {
                matrix3[i] += matrix1[i][j] * matrix2[j];
            }
        }

        return matrix3;
    }

    static mult1D2D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix1D(matrix1) && Matrix3D.checkMatrix2D(matrix2))) {
            return null;
        }

        const matrix3 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix3[i] = 0.0;
            for (let j = 0; j < 3; j++) {
                matrix3[i] += matrix1[j] * matrix2[j][i];
            }
        }

        return matrix3;
    }

    static mult2D2D(matrix1, matrix2) {
        if (!(Matrix3D.checkMatrix2D(matrix1) && Matrix3D.checkMatrix2D(matrix2))) {
            return null;
        }

        const matrix3 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix3[i][j] = 0.0;
                for (let k = 0; k < 3; k++) {
                    matrix3[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }

        return matrix3;
    }

    static multAlpha1D(alpha, matrix) {
        if (!Matrix3D.checkMatrix1D(matrix)) {
            return null;
        }

        const matrix2 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix2[i] = alpha * matrix[i];
        }

        return matrix2;
    }

    static mult1DAlpha(matrix, alpha) {
        return multAlpha1D(alpha, matrix);
    }

    static multAlpha2D(alpha, matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return null;
        }

        const matrix2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix2[i][j] = alpha * matrix[i][j];
            }
        }

        return matrix2;
    }

    static mult2DAlpha(matrix, alpha) {
        return multAlpha2D(alpha, matrix);
    }

    static determinant(matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return 0.0;
        }

        let value = 0.0;
        value += matrix[0][0] * matrix[1][1] * matrix[2][2];
        value += matrix[0][1] * matrix[1][2] * matrix[2][0];
        value += matrix[0][2] * matrix[1][0] * matrix[2][1];
        value -= matrix[0][2] * matrix[1][1] * matrix[2][0];
        value -= matrix[0][0] * matrix[1][2] * matrix[2][1];
        value -= matrix[0][1] * matrix[1][0] * matrix[2][2];

        return value;
    }

    static norm2(matrix) {
        if (!Matrix3D.checkMatrix1D(matrix)) {
            return 0.0;
        }

        let value = 0.0;
        value += matrix[0] * matrix[0];
        value += matrix[1] * matrix[1];
        value += matrix[2] * matrix[2];

        return value;
    }

    static norm(matrix) {
        return Math.sqrt(Matrix3D.norm2(matrix));
    }

    static max1D(matrix) {
        if (!Matrix3D.checkMatrix1D(matrix)) {
            return 0.0;
        }

        return Math.max(Math.max(matrix[0], matrix[1]), matrix[2]);
    }

    static max2D(matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return 0.0;
        }

        let value = matrix[0][0];
        for (let i = 0; i < 3; i++) {
            value = Math.max(value, Math.max(Math.max(matrix[i][0], matrix[i][1]), matrix[i][2]));
        }

        return value;
    }

    static inverse(matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return null;
        }

        const det = Matrix3D.determinant(matrix);
        if (Math.abs(det) < Matrix3D.MIN_DET) {
            return null;
        }

        const matrix2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        matrix2[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
        matrix2[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det;
        matrix2[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;
        matrix2[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det;
        matrix2[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
        matrix2[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det;
        matrix2[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
        matrix2[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det;
        matrix2[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;

        return matrix2;
    }

    static copy2D(matrix) {
        if (!Matrix3D.checkMatrix2D(matrix)) {
            return null;
        }

        const matrix2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                matrix2[i][j] = matrix[i][j];
            }
        }

        return matrix2;
    }

    static copy1D(matrix) {
        if (!Matrix3D.checkMatrix1D(matrix)) {
            return null;
        }

        const matrix2 = [0.0, 0.0, 0.0];
        for (let i = 0; i < 3; i++) {
            matrix2[i] = matrix[i];
        }

        return matrix2;
    }

    static equals2D2D(matrix1, matrix2) {
        return Matrix3D.equals2D2DThreshold(matrix1, matrix2, Matrix3D.EQUALS_THRESHOLD);
    }

    static equals2D2DThreshold(matrix1, matrix2, threshold) {
        if (!(Matrix3D.checkMatrix2D(matrix1) && Matrix3D.checkMatrix2D(matrix2))) {
            return false;
        }

        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                if (Math.abs(matrix1[i][j] - matrix2[i][j]) > threshold) {
                    return false;
                }
            }
        }

        return true;
    }

    static equals1D1D(matrix1, matrix2) {
        return Matrix3D.equals(matrix1, matrix2, EQUALS_THRESHOLD);
    }

    static equals1D1DThreshold(matrix1, matrix2, threshold) {
        if (!(Matrix3D.checkMatrix1D(matrix1) && Matrix3D.checkMatrix1D(matrix2))) {
            return false;
        }

        for (let i = 0; i < 3; i++) {
            if (Math.abs(matrix1[i] - matrix2[i]) > threshold) {
                return false;
            }
        }

        return true;
    }
}
