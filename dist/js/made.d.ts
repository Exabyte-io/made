import { Basis } from "./basis/basis";
import { AtomicConstraints } from "./constraints/constraints";
import { Lattice } from "./lattice/lattice";
import { ReciprocalLattice } from "./lattice/reciprocal/lattice_reciprocal";
import { MaterialMixin } from "./material";
export declare const Made: {
    coefficients: {
        EV_TO_RY: number;
        BOHR_TO_ANGSTROM: number;
        ANGSTROM_TO_BOHR: number;
        EV_A_TO_RY_BOHR: number;
    };
    tolerance: {
        length: number;
        lengthAngstrom: number;
        pointsDistance: number;
    };
    units: {
        bohr: string;
        angstrom: string;
        degree: string;
        radian: string;
        alat: string;
    };
    ATOMIC_COORD_UNITS: {
        crystal: string;
        cartesian: string;
    };
    math: {
        PI: number;
        trunc: (x: number) => number;
        product: (v1: number[], v2: number[]) => number;
        vlen: (v: number[]) => number;
        angle: (a: number[], b: number[], unit: string) => number;
        angleUpTo90: (a: number[], b: number[], unit: string) => number;
        vDist: (v1: number[], v2: number[]) => number | undefined;
        vEqualWithTolerance: (vec1: number[], vec2: number[], tolerance?: number | undefined) => boolean;
        mod: (num: number, tolerance?: number | undefined) => number;
        isBetweenZeroInclusiveAndOne: (number: number, tolerance?: number | undefined) => boolean;
        cartesianProduct: (...arg: number[][]) => number[][];
        almostEqual: (a: number, b: number, tolerance?: number | undefined) => boolean;
        combinations: (a: number, b: number, c: number) => number[][];
        combinationsFromIntervals: (arrA: number[], arrB: number[], arrC: number[]) => number[][];
        calculateSegmentsBetweenPoints3D: (point1: (string | number)[], point2: (string | number)[], n: string | number) => number[][];
        roundToZero: (n: number) => number;
        precise: (x: number, n?: number | undefined) => number;
        roundValueToNDecimals: (value: number, decimals?: number | undefined) => number;
        numberToPrecision: typeof import("@mat3ra/code/dist/js/math").numberToPrecision;
        roundCustom: (value: number, decimals?: number | undefined, method?: import("@mat3ra/code/dist/js/math").RoundingMethodEnum | undefined) => number;
        RoundingMethod: typeof import("@mat3ra/code/dist/js/math").RoundingMethodEnum;
        roundArrayOrNumber: (value: unknown, decimals?: number | undefined, method?: import("@mat3ra/code/dist/js/math").RoundingMethodEnum | undefined) => unknown;
        e: number;
        pi: number;
        i: number;
        Infinity: number;
        LN2: number;
        LN10: number;
        LOG2E: number;
        LOG10E: number;
        NaN: number;
        null: number;
        phi: number;
        SQRT1_2: number;
        SQRT2: number;
        tau: number;
        uninitialized: any;
        version: string;
        expression: import("mathjs").MathNode;
        json: import("mathjs").MathJsJson;
        config: (options: import("mathjs").ConfigOptions) => import("mathjs").ConfigOptions;
        typed: (name: string, signatures: Record<string, (...args: any[]) => any>) => (...args: any[]) => any;
        bignumber(x?: string | number | boolean | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | null | undefined): import("mathjs").BigNumber;
        boolean(x: string | number | boolean | import("mathjs").MathArray | import("mathjs").Matrix | null): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        chain(value?: any): import("mathjs").MathJsChain;
        complex(arg?: string | import("mathjs").Complex | import("mathjs").PolarCoordinates | undefined): import("mathjs").Complex;
        complex(arg?: import("mathjs").MathArray | import("mathjs").Matrix | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        complex(re: number, im: number): import("mathjs").Complex;
        createUnit(name: string, definition?: string | import("mathjs").UnitDefinition | undefined, options?: import("mathjs").CreateUnitOptions | undefined): import("mathjs").Unit;
        createUnit(units: Record<string, string | import("mathjs").UnitDefinition>, options?: import("mathjs").CreateUnitOptions | undefined): import("mathjs").Unit;
        fraction(args: import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Fraction): import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Fraction;
        fraction(numerator: string | number | import("mathjs").MathArray | import("mathjs").Matrix, denominator?: string | number | import("mathjs").MathArray | import("mathjs").Matrix | undefined): import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Fraction;
        index(...ranges: any[]): import("mathjs").Index;
        matrix(format?: "sparse" | "dense" | undefined): import("mathjs").Matrix;
        matrix(data: import("mathjs").MathArray | import("mathjs").Matrix, format?: "sparse" | "dense" | undefined, dataType?: string | undefined): import("mathjs").Matrix;
        number(value?: string | number | boolean | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Unit | null | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix;
        number(unit: import("mathjs").Unit, valuelessUnit: string | import("mathjs").Unit): number;
        sparse(data?: import("mathjs").MathArray | import("mathjs").Matrix | undefined, dataType?: string | undefined): import("mathjs").Matrix;
        splitUnit(unit: import("mathjs").Unit, parts: import("mathjs").Unit[]): import("mathjs").Unit[];
        string(value: import("mathjs").MathType | null): string | import("mathjs").MathArray | import("mathjs").Matrix;
        unit(unit: string): import("mathjs").Unit;
        unit(value: number | import("mathjs").MathArray | import("mathjs").Matrix, unit: string): import("mathjs").Unit;
        compile(expr: import("mathjs").MathExpression): import("mathjs").EvalFunction;
        compile(exprs: import("mathjs").MathExpression[]): import("mathjs").EvalFunction[];
        eval(expr: import("mathjs").MathExpression | import("mathjs").MathExpression[], scope?: object | undefined): any;
        help(search: () => any): import("mathjs").Help;
        parse(expr: import("mathjs").MathExpression, options?: any): import("mathjs").MathNode;
        parse(exprs: import("mathjs").MathExpression[], options?: any): import("mathjs").MathNode[];
        parser(): import("mathjs").Parser;
        derivative(expr: string | import("mathjs").MathNode, variable: string | import("mathjs").MathNode, options?: {
            simplify: boolean;
        } | undefined): import("mathjs").MathNode;
        lsolve(L: import("mathjs").MathArray | import("mathjs").Matrix, b: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        lup(A?: import("mathjs").MathArray | import("mathjs").Matrix | undefined): {
            L: import("mathjs").MathArray | import("mathjs").Matrix;
            U: import("mathjs").MathArray | import("mathjs").Matrix;
            P: number[];
        };
        lusolve(A: number | import("mathjs").MathArray | import("mathjs").Matrix, b: import("mathjs").MathArray | import("mathjs").Matrix, order?: number | undefined, threshold?: number | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        qr(A: import("mathjs").MathArray | import("mathjs").Matrix): {
            Q: import("mathjs").MathArray | import("mathjs").Matrix;
            R: import("mathjs").MathArray | import("mathjs").Matrix;
        };
        rationalize(expr: string | import("mathjs").MathNode, optional?: boolean | object | undefined, detailed?: true | undefined): {
            expression: string | import("mathjs").MathNode;
            variables: string[];
            coefficients: import("mathjs").MathType[];
        };
        rationalize(expr: string | import("mathjs").MathNode, optional?: boolean | object | undefined, detailed?: false | undefined): import("mathjs").MathNode;
        simplify(expr: string | import("mathjs").MathNode, rules?: (string | {
            l: string;
            r: string;
        } | ((node: import("mathjs").MathNode) => import("mathjs").MathNode))[] | undefined, scope?: object | undefined): import("mathjs").MathNode;
        slu(A: import("mathjs").Matrix, order: number, threshold: number): object;
        usolve(U: import("mathjs").MathArray | import("mathjs").Matrix, b: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        abs(x: number): number;
        abs(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        abs(x: import("mathjs").Fraction): import("mathjs").Fraction;
        abs(x: import("mathjs").Complex): import("mathjs").Complex;
        abs(x: import("mathjs").MathArray): import("mathjs").MathArray;
        abs(x: import("mathjs").Matrix): import("mathjs").Matrix;
        abs(x: import("mathjs").Unit): import("mathjs").Unit;
        add(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        cbrt(x: number, allRoots?: boolean | undefined): number;
        cbrt(x: import("mathjs").BigNumber, allRoots?: boolean | undefined): import("mathjs").BigNumber;
        cbrt(x: import("mathjs").Fraction, allRoots?: boolean | undefined): import("mathjs").Fraction;
        cbrt(x: import("mathjs").Complex, allRoots?: boolean | undefined): import("mathjs").Complex;
        cbrt(x: import("mathjs").MathArray, allRoots?: boolean | undefined): import("mathjs").MathArray;
        cbrt(x: import("mathjs").Matrix, allRoots?: boolean | undefined): import("mathjs").Matrix;
        cbrt(x: import("mathjs").Unit, allRoots?: boolean | undefined): import("mathjs").Unit;
        ceil(x: number): number;
        ceil(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        ceil(x: import("mathjs").Fraction): import("mathjs").Fraction;
        ceil(x: import("mathjs").Complex): import("mathjs").Complex;
        ceil(x: import("mathjs").MathArray): import("mathjs").MathArray;
        ceil(x: import("mathjs").Matrix): import("mathjs").Matrix;
        ceil(x: import("mathjs").Unit): import("mathjs").Unit;
        cube(x: number): number;
        cube(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        cube(x: import("mathjs").Fraction): import("mathjs").Fraction;
        cube(x: import("mathjs").Complex): import("mathjs").Complex;
        cube(x: import("mathjs").MathArray): import("mathjs").MathArray;
        cube(x: import("mathjs").Matrix): import("mathjs").Matrix;
        cube(x: import("mathjs").Unit): import("mathjs").Unit;
        divide(x: import("mathjs").Unit, y: import("mathjs").Unit): import("mathjs").Unit;
        divide(x: number, y: number): number;
        divide(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        dotDivide(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        dotMultiply(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        dotPow(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        exp(x: number): number;
        exp(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        exp(x: import("mathjs").Complex): import("mathjs").Complex;
        exp(x: import("mathjs").MathArray): import("mathjs").MathArray;
        exp(x: import("mathjs").Matrix): import("mathjs").Matrix;
        expm1(x: number): number;
        expm1(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        expm1(x: import("mathjs").Complex): import("mathjs").Complex;
        expm1(x: import("mathjs").MathArray): import("mathjs").MathArray;
        expm1(x: import("mathjs").Matrix): import("mathjs").Matrix;
        fix(x: number): number;
        fix(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        fix(x: import("mathjs").Fraction): import("mathjs").Fraction;
        fix(x: import("mathjs").Complex): import("mathjs").Complex;
        fix(x: import("mathjs").MathArray): import("mathjs").MathArray;
        fix(x: import("mathjs").Matrix): import("mathjs").Matrix;
        floor(x: number): number;
        floor(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        floor(x: import("mathjs").Fraction): import("mathjs").Fraction;
        floor(x: import("mathjs").Complex): import("mathjs").Complex;
        floor(x: import("mathjs").MathArray): import("mathjs").MathArray;
        floor(x: import("mathjs").Matrix): import("mathjs").Matrix;
        gcd(...args: number[]): number;
        gcd(...args: import("mathjs").BigNumber[]): import("mathjs").BigNumber;
        gcd(...args: import("mathjs").Fraction[]): import("mathjs").Fraction;
        gcd(...args: import("mathjs").MathArray[]): import("mathjs").MathArray;
        gcd(...args: import("mathjs").Matrix[]): import("mathjs").Matrix;
        hypot(...args: number[]): number;
        hypot(...args: import("mathjs").BigNumber[]): import("mathjs").BigNumber;
        lcm(a: number, b: number): number;
        lcm(a: import("mathjs").BigNumber, b: import("mathjs").BigNumber): import("mathjs").BigNumber;
        lcm(a: import("mathjs").MathArray, b: import("mathjs").MathArray): import("mathjs").MathArray;
        lcm(a: import("mathjs").Matrix, b: import("mathjs").Matrix): import("mathjs").Matrix;
        log(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex;
        log10(x: number): number;
        log10(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        log10(x: import("mathjs").Complex): import("mathjs").Complex;
        log10(x: import("mathjs").MathArray): import("mathjs").MathArray;
        log10(x: import("mathjs").Matrix): import("mathjs").Matrix;
        log1p(x: number, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): number;
        log1p(x: import("mathjs").BigNumber, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): import("mathjs").BigNumber;
        log1p(x: import("mathjs").Complex, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): import("mathjs").Complex;
        log1p(x: import("mathjs").MathArray, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): import("mathjs").MathArray;
        log1p(x: import("mathjs").Matrix, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): import("mathjs").Matrix;
        log2(x: number): number;
        log2(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        log2(x: import("mathjs").Complex): import("mathjs").Complex;
        log2(x: import("mathjs").MathArray): import("mathjs").MathArray;
        log2(x: import("mathjs").Matrix): import("mathjs").Matrix;
        multiply(x: import("mathjs").MathArray | import("mathjs").Matrix, y: import("mathjs").MathType): import("mathjs").MathArray | import("mathjs").Matrix;
        multiply(x: import("mathjs").Unit, y: import("mathjs").Unit): import("mathjs").Unit;
        multiply(x: number, y: number): number;
        multiply(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        norm(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex, p?: string | number | import("mathjs").BigNumber | undefined): number | import("mathjs").BigNumber;
        nthRoot(a: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex, root?: number | import("mathjs").BigNumber | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Complex;
        pow(x: import("mathjs").MathType, y: number | import("mathjs").BigNumber | import("mathjs").Complex): import("mathjs").MathType;
        round(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex, n?: number | import("mathjs").MathArray | import("mathjs").BigNumber | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex;
        sign(x: number): number;
        sign(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        sign(x: import("mathjs").Fraction): import("mathjs").Fraction;
        sign(x: import("mathjs").Complex): import("mathjs").Complex;
        sign(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sign(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sign(x: import("mathjs").Unit): import("mathjs").Unit;
        sqrt(x: number): number;
        sqrt(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        sqrt(x: import("mathjs").Complex): import("mathjs").Complex;
        sqrt(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sqrt(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sqrt(x: import("mathjs").Unit): import("mathjs").Unit;
        square(x: number): number;
        square(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        square(x: import("mathjs").Fraction): import("mathjs").Fraction;
        square(x: import("mathjs").Complex): import("mathjs").Complex;
        square(x: import("mathjs").MathArray): import("mathjs").MathArray;
        square(x: import("mathjs").Matrix): import("mathjs").Matrix;
        square(x: import("mathjs").Unit): import("mathjs").Unit;
        subtract(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
        unaryMinus(x: number): number;
        unaryMinus(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        unaryMinus(x: import("mathjs").Fraction): import("mathjs").Fraction;
        unaryMinus(x: import("mathjs").Complex): import("mathjs").Complex;
        unaryMinus(x: import("mathjs").MathArray): import("mathjs").MathArray;
        unaryMinus(x: import("mathjs").Matrix): import("mathjs").Matrix;
        unaryMinus(x: import("mathjs").Unit): import("mathjs").Unit;
        unaryPlus(x: number): number;
        unaryPlus(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        unaryPlus(x: import("mathjs").Fraction): import("mathjs").Fraction;
        unaryPlus(x: string): string;
        unaryPlus(x: import("mathjs").Complex): import("mathjs").Complex;
        unaryPlus(x: import("mathjs").MathArray): import("mathjs").MathArray;
        unaryPlus(x: import("mathjs").Matrix): import("mathjs").Matrix;
        unaryPlus(x: import("mathjs").Unit): import("mathjs").Unit;
        xgcd(a: number | import("mathjs").BigNumber, b: number | import("mathjs").BigNumber): import("mathjs").MathArray;
        bitAnd(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber, y: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        bitNot(x: number): number;
        bitNot(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        bitNot(x: import("mathjs").MathArray): import("mathjs").MathArray;
        bitNot(x: import("mathjs").Matrix): import("mathjs").Matrix;
        bitOr(x: number, y: number): number;
        bitOr(x: import("mathjs").BigNumber, y: import("mathjs").BigNumber): import("mathjs").BigNumber;
        bitOr(x: import("mathjs").MathArray, y: import("mathjs").MathArray): import("mathjs").MathArray;
        bitOr(x: import("mathjs").Matrix, y: import("mathjs").Matrix): import("mathjs").Matrix;
        bitXor(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber, y: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        leftShift(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber, y: number | import("mathjs").BigNumber): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        rightArithShift(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber, y: number | import("mathjs").BigNumber): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        rightLogShift(x: number | import("mathjs").MathArray | import("mathjs").Matrix, y: number): number | import("mathjs").MathArray | import("mathjs").Matrix;
        bellNumbers(n: number): number;
        bellNumbers(n: import("mathjs").BigNumber): import("mathjs").BigNumber;
        catalan(n: number): number;
        catalan(n: import("mathjs").BigNumber): import("mathjs").BigNumber;
        composition(n: number | import("mathjs").BigNumber, k: number | import("mathjs").BigNumber): number | import("mathjs").BigNumber;
        stirlingS2(n: number | import("mathjs").BigNumber, k: number | import("mathjs").BigNumber): number | import("mathjs").BigNumber;
        arg(x: number | import("mathjs").Complex): number;
        arg(x: import("mathjs").BigNumber | import("mathjs").Complex): import("mathjs").BigNumber;
        arg(x: import("mathjs").MathArray): import("mathjs").MathArray;
        arg(x: import("mathjs").Matrix): import("mathjs").Matrix;
        conj(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex;
        im(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        re(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        distance(x: object | import("mathjs").MathArray | import("mathjs").Matrix, y: object | import("mathjs").MathArray | import("mathjs").Matrix): number | import("mathjs").BigNumber;
        intersect(w: import("mathjs").MathArray | import("mathjs").Matrix, x: import("mathjs").MathArray | import("mathjs").Matrix, y: import("mathjs").MathArray | import("mathjs").Matrix, z: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray;
        and(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit, y: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        not(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        or(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit, y: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        xor(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit, y: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        concat(...args: (import("mathjs").MathArray | import("mathjs").Matrix)[]): import("mathjs").MathArray | import("mathjs").Matrix;
        cross(x: import("mathjs").MathArray | import("mathjs").Matrix, y: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        det(x: import("mathjs").MathArray | import("mathjs").Matrix): number;
        diag(X: import("mathjs").MathArray | import("mathjs").Matrix, format?: string | undefined): import("mathjs").Matrix;
        diag(X: import("mathjs").MathArray | import("mathjs").Matrix, k: number | import("mathjs").BigNumber, format?: string | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        dot(x: import("mathjs").MathArray | import("mathjs").Matrix, y: import("mathjs").MathArray | import("mathjs").Matrix): number;
        expm(x: import("mathjs").Matrix): import("mathjs").Matrix;
        identity(size: number | import("mathjs").MathArray | import("mathjs").Matrix, format?: string | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix;
        identity(m: number, n: number, format?: string | undefined): number | import("mathjs").MathArray | import("mathjs").Matrix;
        filter(x: import("mathjs").MathArray | import("mathjs").Matrix | string[], test: RegExp | ((value: any, index: any, matrix: import("mathjs").MathArray | import("mathjs").Matrix | string[]) => boolean)): import("mathjs").MathArray | import("mathjs").Matrix;
        flatten(x: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        forEach(x: import("mathjs").MathArray | import("mathjs").Matrix, callback: (value: any, index: any, matrix: import("mathjs").MathArray | import("mathjs").Matrix) => void): void;
        inv(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Complex): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Complex;
        kron(x: import("mathjs").MathArray | import("mathjs").Matrix, y: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").Matrix;
        map(x: import("mathjs").MathArray | import("mathjs").Matrix, callback: (value: any, index: any, matrix: import("mathjs").MathArray | import("mathjs").Matrix) => string | import("mathjs").MathType): import("mathjs").MathArray | import("mathjs").Matrix;
        ones(size: number | number[], format?: string | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        ones(m: number, n: number, format?: string | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        partitionSelect(x: import("mathjs").MathArray | import("mathjs").Matrix, k: number, compare?: "asc" | "desc" | ((a: any, b: any) => number) | undefined): any;
        range(str: string, includeEnd?: boolean | undefined): import("mathjs").Matrix;
        range(start: number | import("mathjs").BigNumber, end: number | import("mathjs").BigNumber, includeEnd?: boolean | undefined): import("mathjs").Matrix;
        range(start: number | import("mathjs").BigNumber, end: number | import("mathjs").BigNumber, step: number | import("mathjs").BigNumber, includeEnd?: boolean | undefined): import("mathjs").Matrix;
        reshape(x: import("mathjs").MathArray | import("mathjs").Matrix, sizes: number[]): import("mathjs").MathArray | import("mathjs").Matrix;
        resize(x: import("mathjs").MathArray | import("mathjs").Matrix, size: import("mathjs").MathArray | import("mathjs").Matrix, defaultValue?: string | number | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        size(x: string | number | boolean | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Complex | import("mathjs").Unit): import("mathjs").MathArray | import("mathjs").Matrix;
        sort(x: import("mathjs").MathArray | import("mathjs").Matrix, compare: "asc" | "desc" | ((a: any, b: any) => number) | "natural"): import("mathjs").MathArray | import("mathjs").Matrix;
        sqrtm(A: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        squeeze(x: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        subset(value: string | import("mathjs").MathArray | import("mathjs").Matrix, index: import("mathjs").Index, replacement?: any, defaultValue?: any): string | import("mathjs").MathArray | import("mathjs").Matrix;
        trace(x: import("mathjs").MathArray | import("mathjs").Matrix): number;
        transpose(x: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        zeros(size: number | number[], format?: string | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        zeros(m: number, n: number, format?: string | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        factorial(n: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber;
        gamma(n: number | import("mathjs").MathArray | import("mathjs").Matrix): number | import("mathjs").MathArray | import("mathjs").Matrix;
        kldivergence(q: import("mathjs").MathArray | import("mathjs").Matrix, p: import("mathjs").MathArray | import("mathjs").Matrix): number;
        multinomial(a: number[] | import("mathjs").BigNumber[]): number | import("mathjs").BigNumber;
        permutations(n: number | import("mathjs").BigNumber, k?: number | import("mathjs").BigNumber | undefined): number | import("mathjs").BigNumber;
        pickRandom(array: number[], number?: number | undefined, weights?: number[] | undefined): number;
        random(min?: number | undefined, max?: number | undefined): number;
        random(size: import("mathjs").MathArray | import("mathjs").Matrix, min?: number | undefined, max?: number | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        randomInt(min: number, max?: number | undefined): number;
        randomInt(size: import("mathjs").MathArray | import("mathjs").Matrix, min?: number | undefined, max?: number | undefined): import("mathjs").MathArray | import("mathjs").Matrix;
        compare(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction;
        compareNatural(x: any, y: any): number;
        compareText(x: string | import("mathjs").MathArray | import("mathjs").Matrix, y: string | import("mathjs").MathArray | import("mathjs").Matrix): number | import("mathjs").MathArray | import("mathjs").Matrix;
        deepEqual(x: import("mathjs").MathType, y: import("mathjs").MathType): number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex | import("mathjs").Unit;
        equal(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        equalText(x: string | import("mathjs").MathArray | import("mathjs").Matrix, y: string | import("mathjs").MathArray | import("mathjs").Matrix): number | import("mathjs").MathArray | import("mathjs").Matrix;
        larger(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        largerEq(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        smaller(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        smallerEq(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        unequal(x: string | import("mathjs").MathType, y: string | import("mathjs").MathType): boolean | import("mathjs").MathArray | import("mathjs").Matrix;
        setCartesian(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setDifference(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setDistinct(a: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setIntersect(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setIsSubset(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): boolean;
        setMultiplicity(e: number | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex, a: import("mathjs").MathArray | import("mathjs").Matrix): number;
        setPowerset(a: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setSize(a: import("mathjs").MathArray | import("mathjs").Matrix): number;
        setSymDifference(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        setUnion(a1: import("mathjs").MathArray | import("mathjs").Matrix, a2: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        erf(x: number | import("mathjs").MathArray | import("mathjs").Matrix): number | import("mathjs").MathArray | import("mathjs").Matrix;
        mad(array: import("mathjs").MathArray | import("mathjs").Matrix): any;
        max(...args: import("mathjs").MathType[]): any;
        max(A: import("mathjs").MathArray | import("mathjs").Matrix, dim?: number | undefined): any;
        mean(...args: import("mathjs").MathType[]): any;
        mean(A: import("mathjs").MathArray | import("mathjs").Matrix, dim?: number | undefined): any;
        median(...args: import("mathjs").MathType[]): any;
        min(...args: import("mathjs").MathType[]): any;
        min(A: import("mathjs").MathArray | import("mathjs").Matrix, dim?: number | undefined): any;
        mode(...args: import("mathjs").MathType[]): any;
        prod(...args: import("mathjs").MathType[]): any;
        quantileSeq(A: import("mathjs").MathArray | import("mathjs").Matrix, prob: number | import("mathjs").MathArray | import("mathjs").BigNumber, sorted?: boolean | undefined): number | import("mathjs").MathArray | import("mathjs").BigNumber | import("mathjs").Unit;
        std(array: import("mathjs").MathArray | import("mathjs").Matrix, normalization?: "unbiased" | "uncorrected" | "biased" | undefined): number;
        sum(...args: (number | import("mathjs").BigNumber | import("mathjs").Fraction)[]): any;
        sum(array: import("mathjs").MathArray | import("mathjs").Matrix): any;
        var(...args: (number | import("mathjs").BigNumber | import("mathjs").Fraction)[]): any;
        var(array: import("mathjs").MathArray | import("mathjs").Matrix, normalization?: "unbiased" | "uncorrected" | "biased" | undefined): any;
        format(value: any, options?: number | import("mathjs").FormatOptions | ((item: any) => string) | undefined, callback?: ((value: any) => string) | undefined): string;
        print(template: string, values: any, precision?: number | undefined, options?: number | object | undefined): void;
        acos(x: number): number;
        acos(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acos(x: import("mathjs").Complex): import("mathjs").Complex;
        acos(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acos(x: import("mathjs").Matrix): import("mathjs").Matrix;
        acosh(x: number): number;
        acosh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acosh(x: import("mathjs").Complex): import("mathjs").Complex;
        acosh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acosh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        acot(x: number): number;
        acot(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acot(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acot(x: import("mathjs").Matrix): import("mathjs").Matrix;
        acoth(x: number): number;
        acoth(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acoth(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acoth(x: import("mathjs").Matrix): import("mathjs").Matrix;
        acsc(x: number): number;
        acsc(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acsc(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acsc(x: import("mathjs").Matrix): import("mathjs").Matrix;
        acsch(x: number): number;
        acsch(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        acsch(x: import("mathjs").MathArray): import("mathjs").MathArray;
        acsch(x: import("mathjs").Matrix): import("mathjs").Matrix;
        asec(x: number): number;
        asec(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        asec(x: import("mathjs").MathArray): import("mathjs").MathArray;
        asec(x: import("mathjs").Matrix): import("mathjs").Matrix;
        asech(x: number): number;
        asech(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        asech(x: import("mathjs").MathArray): import("mathjs").MathArray;
        asech(x: import("mathjs").Matrix): import("mathjs").Matrix;
        asin(x: number): number;
        asin(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        asin(x: import("mathjs").Complex): import("mathjs").Complex;
        asin(x: import("mathjs").MathArray): import("mathjs").MathArray;
        asin(x: import("mathjs").Matrix): import("mathjs").Matrix;
        asinh(x: number): number;
        asinh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        asinh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        asinh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        atan(x: number): number;
        atan(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        atan(x: import("mathjs").MathArray): import("mathjs").MathArray;
        atan(x: import("mathjs").Matrix): import("mathjs").Matrix;
        atan2(y: number, x: number): number;
        atan2(y: import("mathjs").MathArray | import("mathjs").Matrix, x: import("mathjs").MathArray | import("mathjs").Matrix): import("mathjs").MathArray | import("mathjs").Matrix;
        atanh(x: number): number;
        atanh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        atanh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        atanh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        cos(x: number | import("mathjs").Unit): number;
        cos(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        cos(x: import("mathjs").Complex): import("mathjs").Complex;
        cos(x: import("mathjs").MathArray): import("mathjs").MathArray;
        cos(x: import("mathjs").Matrix): import("mathjs").Matrix;
        cosh(x: number | import("mathjs").Unit): number;
        cosh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        cosh(x: import("mathjs").Complex): import("mathjs").Complex;
        cosh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        cosh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        cot(x: number | import("mathjs").Unit): number;
        cot(x: import("mathjs").Complex): import("mathjs").Complex;
        cot(x: import("mathjs").MathArray): import("mathjs").MathArray;
        cot(x: import("mathjs").Matrix): import("mathjs").Matrix;
        coth(x: number | import("mathjs").Unit): number;
        coth(x: import("mathjs").Complex): import("mathjs").Complex;
        coth(x: import("mathjs").MathArray): import("mathjs").MathArray;
        coth(x: import("mathjs").Matrix): import("mathjs").Matrix;
        csc(x: number | import("mathjs").Unit): number;
        csc(x: import("mathjs").Complex): import("mathjs").Complex;
        csc(x: import("mathjs").MathArray): import("mathjs").MathArray;
        csc(x: import("mathjs").Matrix): import("mathjs").Matrix;
        csch(x: number | import("mathjs").Unit): number;
        csch(x: import("mathjs").Complex): import("mathjs").Complex;
        csch(x: import("mathjs").MathArray): import("mathjs").MathArray;
        csch(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sec(x: number | import("mathjs").Unit): number;
        sec(x: import("mathjs").Complex): import("mathjs").Complex;
        sec(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sec(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sech(x: number | import("mathjs").Unit): number;
        sech(x: import("mathjs").Complex): import("mathjs").Complex;
        sech(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sech(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sin(x: number | import("mathjs").Unit): number;
        sin(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        sin(x: import("mathjs").Complex): import("mathjs").Complex;
        sin(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sin(x: import("mathjs").Matrix): import("mathjs").Matrix;
        sinh(x: number | import("mathjs").Unit): number;
        sinh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        sinh(x: import("mathjs").Complex): import("mathjs").Complex;
        sinh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        sinh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        tan(x: number | import("mathjs").Unit): number;
        tan(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        tan(x: import("mathjs").Complex): import("mathjs").Complex;
        tan(x: import("mathjs").MathArray): import("mathjs").MathArray;
        tan(x: import("mathjs").Matrix): import("mathjs").Matrix;
        tanh(x: number | import("mathjs").Unit): number;
        tanh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
        tanh(x: import("mathjs").Complex): import("mathjs").Complex;
        tanh(x: import("mathjs").MathArray): import("mathjs").MathArray;
        tanh(x: import("mathjs").Matrix): import("mathjs").Matrix;
        to(x: import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Unit, unit: string | import("mathjs").Unit): import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").Unit;
        clone(x: any): any;
        isInteger(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction): boolean;
        isNaN(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Unit): boolean;
        isNegative(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Unit): boolean;
        isNumeric(x: any): x is number | boolean | import("mathjs").BigNumber | import("mathjs").Fraction;
        isPositive(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Unit): boolean;
        isPrime(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber): boolean;
        isZero(x: number | import("mathjs").MathArray | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex | import("mathjs").Unit): boolean;
        typeof(x: any): string;
        import(object: import("mathjs").ImportObject | import("mathjs").ImportObject[], options: import("mathjs").ImportOptions): void;
    };
    Material: {
        new (...config: any[]): {
            consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
            _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
            prop: {
                <T = undefined>(name: string, defaultValue: T): T;
                <T_1 = undefined>(name: string): T_1 | undefined;
            } & {
                <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                <T_3 = undefined>(name: string): T_3 | undefined;
            } & {
                <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                <T_5 = undefined>(name: string): T_5 | undefined;
            };
            setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
            unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
            setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
            toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
            toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
            validate: (() => void) & (() => void) & (() => void);
            clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
            isValid: (() => boolean) & (() => boolean) & (() => boolean);
            readonly cls: string;
            getClsName: (() => string) & (() => string) & (() => string);
            getAsEntityReference: {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            } & {
                (byIdOnly: true): {
                    _id: string;
                };
                (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
            };
            getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
            id: string;
            _id: string;
            schemaVersion: string;
            systemName: string;
            readonly slug: string;
            readonly isSystemEntity: boolean;
            metadata: object;
            updateMetadata(object: object): void;
            isDefault: boolean;
            name: string;
            setName(name: string): void;
            src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
            isNonPeriodic: boolean;
            readonly formula: string;
            readonly unitCellFormula: string;
            readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
            readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
            readonly uniqueElements: string[];
            lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            readonly Lattice: Lattice;
            hash: string;
            readonly scaledHash: string;
            calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
            getInchiStringForHash(): string;
            getDerivedPropertyByName(name: string): {
                name?: "volume" | undefined;
                units?: "angstrom^3" | undefined;
                value: number;
            } | {
                name?: "density" | undefined;
                units?: "g/cm^3" | undefined;
                value: number;
            } | {
                pointGroupSymbol?: string | undefined;
                spaceGroupSymbol?: string | undefined;
                tolerance?: {
                    units?: "angstrom" | undefined;
                    value: number;
                } | undefined;
                name?: "symmetry" | undefined;
            } | {
                name?: "elemental_ratio" | undefined;
                value: number;
                element?: string | undefined;
            } | {
                name?: "p-norm" | undefined;
                degree?: number | undefined;
                value: number;
            } | {
                name?: "inchi" | undefined;
                value: string;
            } | {
                name?: "inchi_key" | undefined;
                value: string;
            } | undefined;
            getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
            unsetFileProps(): void;
            setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
            updateFormula(): void;
            setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
            toCrystal(): void;
            toCartesian(): void;
            getBasisAsXyz(fractional?: boolean): string;
            getAsQEFormat(): string;
            getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
            getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                name: string;
            } & {
                setName(name: string): void;
            };
            getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
            readonly defaultConfig: {
                name: string;
                basis: {
                    elements: {
                        id: number;
                        value: string;
                    }[];
                    coordinates: {
                        id: number;
                        value: number[];
                    }[];
                    units: string;
                };
                lattice: {
                    type: string;
                    a: number;
                    b: number;
                    c: number;
                    alpha: number;
                    beta: number;
                    gamma: number;
                    units: {
                        length: string;
                        angle: string;
                    };
                };
            };
            constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
        };
    } & (new (...args: any[]) => {
        consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
        addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
        _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        prop<T = undefined>(name: string, defaultValue: T): T;
        prop<T_1 = undefined>(name: string): T_1 | undefined;
        setProp(name: string, value: unknown): void;
        unsetProp(name: string): void;
        setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
        toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        clone(extraContext?: object | undefined): any;
        validate(): void;
        clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        isValid(): boolean;
        readonly cls: string;
        getClsName(): string;
        getAsEntityReference(byIdOnly: true): {
            _id: string;
        };
        getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
        id: string;
        _id: string;
        schemaVersion: string;
        systemName: string;
        readonly slug: string;
        readonly isSystemEntity: boolean;
    }) & (new (...args: any[]) => {
        metadata: object;
        updateMetadata(object: object): void;
        _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
        prop<T_3 = undefined>(name: string): T_3 | undefined;
        setProp(name: string, value: unknown): void;
        unsetProp(name: string): void;
        setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
        toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        clone(extraContext?: object | undefined): any;
        validate(): void;
        clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
        isValid(): boolean;
        readonly cls: string;
        getClsName(): string;
        getAsEntityReference(byIdOnly: true): {
            _id: string;
        };
        getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
        getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
        id: string;
        _id: string;
        schemaVersion: string;
        systemName: string;
        readonly slug: string;
        readonly isSystemEntity: boolean;
    }) & typeof import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & import("@mat3ra/code/dist/js/utils/types").Constructor<{
        isDefault: boolean;
    }> & {
        createDefault(this: import("@mat3ra/code/dist/js/utils/types").Constructor<import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity> & {
            defaultConfig?: object | null | undefined;
        }): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
    } & import("@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin").NamedEntityConstructor;
    MaterialMixin: typeof MaterialMixin;
    defaultMaterialConfig: {
        name: string;
        basis: {
            elements: {
                id: number;
                value: string;
            }[];
            coordinates: {
                id: number;
                value: number[];
            }[];
            units: string;
        };
        lattice: {
            type: string;
            a: number;
            b: number;
            c: number;
            alpha: number;
            beta: number;
            gamma: number;
            units: {
                length: string;
                angle: string;
            };
        };
    };
    Lattice: typeof Lattice;
    nonPeriodicLatticeScalingFactor: number;
    ReciprocalLattice: typeof ReciprocalLattice;
    Basis: typeof Basis;
    AtomicConstraints: typeof AtomicConstraints;
    parsers: {
        xyz: {
            validate: typeof import("./parsers/xyz").validate;
            fromMaterial: (materialOrConfig: import("@mat3ra/esse/dist/js/types").MaterialSchema, fractional?: boolean) => string;
            toBasisConfig: (txt: string, units?: string, cell?: import("./cell/cell").Cell) => import("./basis/constrained_basis").ConstrainedBasisConfig;
            fromBasis: (basisClsInstance: import("./basis/constrained_basis").ConstrainedBasis, coordinatePrintFormat: string) => string;
            CombinatorialBasis: typeof import("./parsers/xyz_combinatorial_basis").CombinatorialBasis;
        };
        poscar: {
            isPoscar: (text: string) => boolean;
            toPoscar: (materialOrConfig: import("./types").MaterialJSON, omitConstraints?: boolean) => string;
            fromPoscar: (fileContent: string) => object;
            atomicConstraintsCharFromBool: (bool: boolean) => string;
            atomsCount: typeof import("./parsers/poscar").atomsCount;
        };
        cif: {
            parseMeta: (txt: string) => import("./parsers/cif").Meta;
        };
        espresso: {
            toEspressoFormat: (materialOrConfig: import("@mat3ra/esse/dist/js/types").MaterialSchema) => string;
        };
        nativeFormatParsers: {
            detectFormat: (text: string) => string;
            convertFromNativeFormat: (text: string) => any;
        };
    };
    tools: {
        surface: {
            generateConfig: (material: {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
                prop: {
                    <T = undefined>(name: string, defaultValue: T): T;
                    <T_1 = undefined>(name: string): T_1 | undefined;
                } & {
                    <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                    <T_3 = undefined>(name: string): T_3 | undefined;
                } & {
                    <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                    <T_5 = undefined>(name: string): T_5 | undefined;
                };
                setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
                unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
                setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
                toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
                toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
                validate: (() => void) & (() => void) & (() => void);
                clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                isValid: (() => boolean) & (() => boolean) & (() => boolean);
                readonly cls: string;
                getClsName: (() => string) & (() => string) & (() => string);
                getAsEntityReference: {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                };
                getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
                metadata: object;
                updateMetadata(object: object): void;
                isDefault: boolean;
                name: string;
                setName(name: string): void;
                src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
                isNonPeriodic: boolean;
                readonly formula: string;
                readonly unitCellFormula: string;
                readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
                readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
                readonly uniqueElements: string[];
                lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
                readonly Lattice: Lattice;
                hash: string;
                readonly scaledHash: string;
                calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
                getInchiStringForHash(): string;
                getDerivedPropertyByName(name: string): {
                    name?: "volume" | undefined;
                    units?: "angstrom^3" | undefined;
                    value: number;
                } | {
                    name?: "density" | undefined;
                    units?: "g/cm^3" | undefined;
                    value: number;
                } | {
                    pointGroupSymbol?: string | undefined;
                    spaceGroupSymbol?: string | undefined;
                    tolerance?: {
                        units?: "angstrom" | undefined;
                        value: number;
                    } | undefined;
                    name?: "symmetry" | undefined;
                } | {
                    name?: "elemental_ratio" | undefined;
                    value: number;
                    element?: string | undefined;
                } | {
                    name?: "p-norm" | undefined;
                    degree?: number | undefined;
                    value: number;
                } | {
                    name?: "inchi" | undefined;
                    value: string;
                } | {
                    name?: "inchi_key" | undefined;
                    value: string;
                } | undefined;
                getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
                unsetFileProps(): void;
                setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
                updateFormula(): void;
                setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
                toCrystal(): void;
                toCartesian(): void;
                getBasisAsXyz(fractional?: boolean): string;
                getAsQEFormat(): string;
                getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
                getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                    name: string;
                } & {
                    setName(name: string): void;
                };
                getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                readonly defaultConfig: {
                    name: string;
                    basis: {
                        elements: {
                            id: number;
                            value: string;
                        }[];
                        coordinates: {
                            id: number;
                            value: number[];
                        }[];
                        units: string;
                    };
                    lattice: {
                        type: string;
                        a: number;
                        b: number;
                        c: number;
                        alpha: number;
                        beta: number;
                        gamma: number;
                        units: {
                            length: string;
                            angle: string;
                        };
                    };
                };
                constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
            } & {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T = undefined>(name: string, defaultValue: T): T;
                prop<T_1 = undefined>(name: string): T_1 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & {
                metadata: object;
                updateMetadata(object: object): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                prop<T_3 = undefined>(name: string): T_3 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                isDefault: boolean;
            } & {
                name: string;
            } & {
                setName(name: string): void;
            }, millerIndices: import("@mat3ra/esse/dist/js/types").Coordinate3DSchema, numberOfLayers?: number, vx?: number, vy?: number) => import("./tools/surface").SlabConfigSchema;
        };
        supercell: {
            generateConfig: (material: import("./types").MaterialInterface, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => {
                name: string;
                basis: import("@mat3ra/esse/dist/js/types").BasisSchema;
                lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
            };
            generateNewBasisWithinSupercell: (basis: Basis | import("./basis/constrained_basis").ConstrainedBasis, cell: import("./cell/cell").Cell, supercell: import("./cell/cell").Cell, supercellMatrix: import("@mat3ra/esse/dist/js/types").Matrix3X3Schema) => Basis;
        };
        material: {
            scaleOneLatticeVector: (material: {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
                prop: {
                    <T = undefined>(name: string, defaultValue: T): T;
                    <T_1 = undefined>(name: string): T_1 | undefined;
                } & {
                    <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                    <T_3 = undefined>(name: string): T_3 | undefined;
                } & {
                    <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                    <T_5 = undefined>(name: string): T_5 | undefined;
                };
                setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
                unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
                setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
                toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
                toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
                validate: (() => void) & (() => void) & (() => void);
                clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                isValid: (() => boolean) & (() => boolean) & (() => boolean);
                readonly cls: string;
                getClsName: (() => string) & (() => string) & (() => string);
                getAsEntityReference: {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                };
                getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
                metadata: object;
                updateMetadata(object: object): void;
                isDefault: boolean;
                name: string;
                setName(name: string): void;
                src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
                isNonPeriodic: boolean;
                readonly formula: string;
                readonly unitCellFormula: string;
                readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
                readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
                readonly uniqueElements: string[];
                lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
                readonly Lattice: Lattice;
                hash: string;
                readonly scaledHash: string;
                calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
                getInchiStringForHash(): string;
                getDerivedPropertyByName(name: string): {
                    name?: "volume" | undefined;
                    units?: "angstrom^3" | undefined;
                    value: number;
                } | {
                    name?: "density" | undefined;
                    units?: "g/cm^3" | undefined;
                    value: number;
                } | {
                    pointGroupSymbol?: string | undefined;
                    spaceGroupSymbol?: string | undefined;
                    tolerance?: {
                        units?: "angstrom" | undefined;
                        value: number;
                    } | undefined;
                    name?: "symmetry" | undefined;
                } | {
                    name?: "elemental_ratio" | undefined;
                    value: number;
                    element?: string | undefined;
                } | {
                    name?: "p-norm" | undefined;
                    degree?: number | undefined;
                    value: number;
                } | {
                    name?: "inchi" | undefined;
                    value: string;
                } | {
                    name?: "inchi_key" | undefined;
                    value: string;
                } | undefined;
                getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
                unsetFileProps(): void;
                setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
                updateFormula(): void;
                setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
                toCrystal(): void;
                toCartesian(): void;
                getBasisAsXyz(fractional?: boolean): string;
                getAsQEFormat(): string;
                getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
                getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                    name: string;
                } & {
                    setName(name: string): void;
                };
                getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                readonly defaultConfig: {
                    name: string;
                    basis: {
                        elements: {
                            id: number;
                            value: string;
                        }[];
                        coordinates: {
                            id: number;
                            value: number[];
                        }[];
                        units: string;
                    };
                    lattice: {
                        type: string;
                        a: number;
                        b: number;
                        c: number;
                        alpha: number;
                        beta: number;
                        gamma: number;
                        units: {
                            length: string;
                            angle: string;
                        };
                    };
                };
                constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
            } & {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T = undefined>(name: string, defaultValue: T): T;
                prop<T_1 = undefined>(name: string): T_1 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & {
                metadata: object;
                updateMetadata(object: object): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                prop<T_3 = undefined>(name: string): T_3 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                isDefault: boolean;
            } & {
                name: string;
            } & {
                setName(name: string): void;
            }, key?: "a" | "b" | "c", factor?: number) => void;
            scaleLatticeToMakeNonPeriodic: (material: {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
                prop: {
                    <T = undefined>(name: string, defaultValue: T): T;
                    <T_1 = undefined>(name: string): T_1 | undefined;
                } & {
                    <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                    <T_3 = undefined>(name: string): T_3 | undefined;
                } & {
                    <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                    <T_5 = undefined>(name: string): T_5 | undefined;
                };
                setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
                unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
                setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
                toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
                toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
                validate: (() => void) & (() => void) & (() => void);
                clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                isValid: (() => boolean) & (() => boolean) & (() => boolean);
                readonly cls: string;
                getClsName: (() => string) & (() => string) & (() => string);
                getAsEntityReference: {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                };
                getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
                metadata: object;
                updateMetadata(object: object): void;
                isDefault: boolean;
                name: string;
                setName(name: string): void;
                src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
                isNonPeriodic: boolean;
                readonly formula: string;
                readonly unitCellFormula: string;
                readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
                readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
                readonly uniqueElements: string[];
                lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
                readonly Lattice: Lattice;
                hash: string;
                readonly scaledHash: string;
                calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
                getInchiStringForHash(): string;
                getDerivedPropertyByName(name: string): {
                    name?: "volume" | undefined;
                    units?: "angstrom^3" | undefined;
                    value: number;
                } | {
                    name?: "density" | undefined;
                    units?: "g/cm^3" | undefined;
                    value: number;
                } | {
                    pointGroupSymbol?: string | undefined;
                    spaceGroupSymbol?: string | undefined;
                    tolerance?: {
                        units?: "angstrom" | undefined;
                        value: number;
                    } | undefined;
                    name?: "symmetry" | undefined;
                } | {
                    name?: "elemental_ratio" | undefined;
                    value: number;
                    element?: string | undefined;
                } | {
                    name?: "p-norm" | undefined;
                    degree?: number | undefined;
                    value: number;
                } | {
                    name?: "inchi" | undefined;
                    value: string;
                } | {
                    name?: "inchi_key" | undefined;
                    value: string;
                } | undefined;
                getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
                unsetFileProps(): void;
                setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
                updateFormula(): void;
                setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
                toCrystal(): void;
                toCartesian(): void;
                getBasisAsXyz(fractional?: boolean): string;
                getAsQEFormat(): string;
                getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
                getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                    name: string;
                } & {
                    setName(name: string): void;
                };
                getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                readonly defaultConfig: {
                    name: string;
                    basis: {
                        elements: {
                            id: number;
                            value: string;
                        }[];
                        coordinates: {
                            id: number;
                            value: number[];
                        }[];
                        units: string;
                    };
                    lattice: {
                        type: string;
                        a: number;
                        b: number;
                        c: number;
                        alpha: number;
                        beta: number;
                        gamma: number;
                        units: {
                            length: string;
                            angle: string;
                        };
                    };
                };
                constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
            } & {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T = undefined>(name: string, defaultValue: T): T;
                prop<T_1 = undefined>(name: string): T_1 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & {
                metadata: object;
                updateMetadata(object: object): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                prop<T_3 = undefined>(name: string): T_3 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                isDefault: boolean;
            } & {
                name: string;
            } & {
                setName(name: string): void;
            }) => void;
            getBasisConfigTranslatedToCenter: (material: {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject & import("./materialMixin").MaterialSchemaJSON;
                prop: {
                    <T = undefined>(name: string, defaultValue: T): T;
                    <T_1 = undefined>(name: string): T_1 | undefined;
                } & {
                    <T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                    <T_3 = undefined>(name: string): T_3 | undefined;
                } & {
                    <T_4 = undefined>(name: string, defaultValue: T_4): T_4;
                    <T_5 = undefined>(name: string): T_5 | undefined;
                };
                setProp: ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void) & ((name: string, value: unknown) => void);
                unsetProp: ((name: string) => void) & ((name: string) => void) & ((name: string) => void);
                setProps: ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any) & ((json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined) => any);
                toJSON: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & (() => import("./types").MaterialJSON);
                toJSONSafe: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                toJSONQuick: ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((exclude?: string[] | undefined) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                clone: ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any) & ((extraContext?: object | undefined) => any);
                validate: (() => void) & (() => void) & (() => void);
                clean: ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject) & ((config: import("@mat3ra/esse/dist/js/esse/types").AnyObject) => import("@mat3ra/esse/dist/js/esse/types").AnyObject);
                isValid: (() => boolean) & (() => boolean) & (() => boolean);
                readonly cls: string;
                getClsName: (() => string) & (() => string) & (() => string);
                getAsEntityReference: {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                } & {
                    (byIdOnly: true): {
                        _id: string;
                    };
                    (byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                };
                getEntityByName: ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity) & ((entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string) => import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity);
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
                metadata: object;
                updateMetadata(object: object): void;
                isDefault: boolean;
                name: string;
                setName(name: string): void;
                src: import("@mat3ra/esse/dist/js/types").FileSourceSchema | undefined;
                isNonPeriodic: boolean;
                readonly formula: string;
                readonly unitCellFormula: string;
                readonly basis: import("./materialMixin").OptionallyConstrainedBasisConfig;
                readonly Basis: import("./basis/constrained_basis").ConstrainedBasis;
                readonly uniqueElements: string[];
                lattice: import("@mat3ra/esse/dist/js/types").LatticeSchema;
                readonly Lattice: Lattice;
                hash: string;
                readonly scaledHash: string;
                calculateHash(salt?: string, isScaled?: boolean, bypassNonPeriodicCheck?: boolean): string;
                getInchiStringForHash(): string;
                getDerivedPropertyByName(name: string): {
                    name?: "volume" | undefined;
                    units?: "angstrom^3" | undefined;
                    value: number;
                } | {
                    name?: "density" | undefined;
                    units?: "g/cm^3" | undefined;
                    value: number;
                } | {
                    pointGroupSymbol?: string | undefined;
                    spaceGroupSymbol?: string | undefined;
                    tolerance?: {
                        units?: "angstrom" | undefined;
                        value: number;
                    } | undefined;
                    name?: "symmetry" | undefined;
                } | {
                    name?: "elemental_ratio" | undefined;
                    value: number;
                    element?: string | undefined;
                } | {
                    name?: "p-norm" | undefined;
                    degree?: number | undefined;
                    value: number;
                } | {
                    name?: "inchi" | undefined;
                    value: string;
                } | {
                    name?: "inchi_key" | undefined;
                    value: string;
                } | undefined;
                getDerivedProperties(): import("@mat3ra/esse/dist/js/types").DerivedPropertiesSchema;
                unsetFileProps(): void;
                setBasis(textOrObject: string | import("./basis/basis").BasisConfig, format?: string | undefined, unitz?: string | undefined): void;
                updateFormula(): void;
                setBasisConstraints(constraints: import("./constraints/constraints").Constraint[]): void;
                toCrystal(): void;
                toCartesian(): void;
                getBasisAsXyz(fractional?: boolean): string;
                getAsQEFormat(): string;
                getAsPOSCAR(ignoreOriginal?: boolean, omitConstraints?: boolean): string;
                getACopyWithConventionalCell(): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                    name: string;
                } & {
                    setName(name: string): void;
                };
                getConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                getBasisConsistencyChecks(): import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                readonly defaultConfig: {
                    name: string;
                    basis: {
                        elements: {
                            id: number;
                            value: string;
                        }[];
                        coordinates: {
                            id: number;
                            value: number[];
                        }[];
                        units: string;
                    };
                    lattice: {
                        type: string;
                        a: number;
                        b: number;
                        c: number;
                        alpha: number;
                        beta: number;
                        gamma: number;
                        units: {
                            length: string;
                            angle: string;
                        };
                    };
                };
                constructMaterialFileSource(fileName: string, fileContent: string, fileExtension: string): import("@mat3ra/esse/dist/js/types").FileSourceSchema;
            } & {
                consistencyChecks: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[];
                addConsistencyChecks(array: import("@mat3ra/esse/dist/js/types").ConsistencyCheck[]): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T = undefined>(name: string, defaultValue: T): T;
                prop<T_1 = undefined>(name: string): T_1 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & {
                metadata: object;
                updateMetadata(object: object): void;
                _json: import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                prop<T_2 = undefined>(name: string, defaultValue: T_2): T_2;
                prop<T_3 = undefined>(name: string): T_3 | undefined;
                setProp(name: string, value: unknown): void;
                unsetProp(name: string): void;
                setProps(json?: import("@mat3ra/esse/dist/js/esse/types").AnyObject | undefined): any;
                toJSON(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONSafe(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                toJSONQuick(exclude?: string[] | undefined): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                clone(extraContext?: object | undefined): any;
                validate(): void;
                clean(config: import("@mat3ra/esse/dist/js/esse/types").AnyObject): import("@mat3ra/esse/dist/js/esse/types").AnyObject;
                isValid(): boolean;
                readonly cls: string;
                getClsName(): string;
                getAsEntityReference(byIdOnly: true): {
                    _id: string;
                };
                getAsEntityReference(byIdOnly?: false | undefined): Required<import("@mat3ra/esse/dist/js/types").EntityReferenceSchema>;
                getEntityByName(entities: import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity[], entity: string, name: string): import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity;
                id: string;
                _id: string;
                schemaVersion: string;
                systemName: string;
                readonly slug: string;
                readonly isSystemEntity: boolean;
            } & import("@mat3ra/code/dist/js/entity/in_memory").InMemoryEntity & {
                isDefault: boolean;
            } & {
                name: string;
            } & {
                setName(name: string): void;
            }) => void;
        };
        basis: {
            repeat: (basis: Basis, repetitions: number[]) => Basis;
            interpolate: (initialBasis: Basis, finalBasis: Basis, numberOfSteps?: number) => Basis[];
        };
    };
    LATTICE_TYPE_CONFIGS: import("./lattice/lattice_types").LatticeTypeConfig[];
    DEFAULT_LATTICE_UNITS: {
        length: {
            angstrom: string;
        };
        angle: {
            degree: string;
        };
    };
};
