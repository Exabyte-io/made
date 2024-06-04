declare const _default: {
    PI: number;
    trunc: (x: number) => number;
    product: (v1: number[], v2: number[]) => import("mathjs").Matrix;
    vlen: (v: number[]) => import("mathjs").Matrix;
    angle: (a: number[], b: number[], unit: string) => number;
    angleUpTo90: (a: number[], b: number[], unit: string) => number;
    vDist: (v1: number[], v2: number[]) => import("mathjs").Matrix | undefined;
    vEqualWithTolerance: (vec1: number[], vec2: number[], tolerance?: number | undefined) => boolean;
    roundToZero: (n: number) => number;
    precise: (x: number, n?: number | undefined) => number;
    mod: (num: number, tolerance?: number | undefined) => number;
    isBetweenZeroInclusiveAndOne: (number: number, tolerance?: number | undefined) => boolean;
    cartesianProduct: (...arg: number[][]) => number[][];
    almostEqual: (a: number, b: number, tolerance?: number | undefined) => boolean;
    combinations: (a: number, b: number, c: number) => number[][];
    combinationsFromIntervals: (arrA: number[], arrB: number[], arrC: number[]) => number[][];
    calculateSegmentsBetweenPoints3D: (point1: (string | number)[], point2: (string | number)[], n: string | number) => number[][];
    roundValueToNDecimals: (value: number, decimals?: number | undefined) => number;
    numberToPrecision: typeof import("@mat3ra/code/dist/js/math").numberToPrecision;
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
    config: (options: any) => void;
    lsolve(L: import("mathjs").Matrix | import("mathjs").MathArray, b: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    lup(A?: import("mathjs").Matrix | import("mathjs").MathArray | undefined): import("mathjs").MathArray;
    lusolve(A: number | import("mathjs").Matrix | import("mathjs").MathArray, b: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    slu(A: import("mathjs").Matrix, order: number, threshold: number): any;
    usolve(U: import("mathjs").Matrix | import("mathjs").MathArray, b: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
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
    kron(x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix;
    lcm(a: number, b: number): number;
    lcm(a: import("mathjs").BigNumber, b: import("mathjs").BigNumber): import("mathjs").BigNumber;
    lcm(a: import("mathjs").MathArray, b: import("mathjs").MathArray): import("mathjs").MathArray;
    lcm(a: import("mathjs").Matrix, b: import("mathjs").Matrix): import("mathjs").Matrix;
    log(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray, base?: number | import("mathjs").BigNumber | import("mathjs").Complex | undefined): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray;
    log10(x: number): number;
    log10(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
    log10(x: import("mathjs").Complex): import("mathjs").Complex;
    log10(x: import("mathjs").MathArray): import("mathjs").MathArray;
    log10(x: import("mathjs").Matrix): import("mathjs").Matrix;
    multiply(x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").MathType): import("mathjs").Matrix;
    multiply(x: import("mathjs").Unit, y: import("mathjs").Unit): import("mathjs").Unit;
    multiply(x: number, y: number): number;
    multiply(x: import("mathjs").MathType, y: import("mathjs").MathType): import("mathjs").MathType;
    norm(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray, p?: string | number | import("mathjs").BigNumber | undefined): number | import("mathjs").BigNumber;
    nthRoot(a: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray, root?: number | import("mathjs").BigNumber | undefined): number | import("mathjs").Matrix | import("mathjs").Complex | import("mathjs").MathArray;
    pow(x: import("mathjs").MathType, y: number | import("mathjs").BigNumber | import("mathjs").Complex): import("mathjs").MathType;
    round(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex | import("mathjs").MathArray, n?: number | import("mathjs").BigNumber | import("mathjs").MathArray | undefined): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex | import("mathjs").MathArray;
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
    bitAnd(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray, y: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    bitNot(x: number): number;
    bitNot(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
    bitNot(x: import("mathjs").MathArray): import("mathjs").MathArray;
    bitNot(x: import("mathjs").Matrix): import("mathjs").Matrix;
    bitOr(x: number): number;
    bitOr(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
    bitOr(x: import("mathjs").MathArray): import("mathjs").MathArray;
    bitOr(x: import("mathjs").Matrix): import("mathjs").Matrix;
    bitXor(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray, y: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    leftShift(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray, y: number | import("mathjs").BigNumber): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    rightArithShift(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray, y: number | import("mathjs").BigNumber): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    rightLogShift(x: number | import("mathjs").Matrix | import("mathjs").MathArray, y: number): number | import("mathjs").Matrix | import("mathjs").MathArray;
    bellNumbers(n: number): number;
    bellNumbers(n: import("mathjs").BigNumber): import("mathjs").BigNumber;
    catalan(n: number): number;
    catalan(n: import("mathjs").BigNumber): import("mathjs").BigNumber;
    composition(n: number | import("mathjs").BigNumber, k: number | import("mathjs").BigNumber): number | import("mathjs").BigNumber;
    stirlingS2(n: number | import("mathjs").BigNumber, k: number | import("mathjs").BigNumber): number | import("mathjs").BigNumber;
    arg(x: number | import("mathjs").Complex): number;
    arg(x: import("mathjs").MathArray): import("mathjs").MathArray;
    arg(x: import("mathjs").Matrix): import("mathjs").Matrix;
    conj(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray;
    im(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    re(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    bignumber(x?: string | number | boolean | import("mathjs").Matrix | import("mathjs").MathArray | undefined): import("mathjs").BigNumber;
    boolean(x: string | number | boolean | import("mathjs").Matrix | import("mathjs").MathArray): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    chain(value?: any): import("mathjs").MathJsChain;
    complex(arg?: string | import("mathjs").Complex | import("mathjs").MathArray | import("mathjs").PolarCoordinates | undefined): import("mathjs").Complex;
    complex(re: number, im: number): import("mathjs").Complex;
    fraction(numerator: string | number | import("mathjs").Matrix | import("mathjs").MathArray, denominator?: string | number | import("mathjs").Matrix | import("mathjs").MathArray | undefined): import("mathjs").Matrix | import("mathjs").Fraction | import("mathjs").MathArray;
    index(...ranges: any[]): import("mathjs").Index;
    matrix(format?: "sparse" | "dense" | undefined): import("mathjs").Matrix;
    matrix(data: import("mathjs").Matrix | import("mathjs").MathArray, format?: "sparse" | "dense" | undefined, dataType?: string | undefined): import("mathjs").Matrix;
    number(value?: string | number | boolean | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Unit | import("mathjs").MathArray | undefined): number | import("mathjs").Matrix | import("mathjs").MathArray;
    number(unit: import("mathjs").Unit, valuelessUnit: string | import("mathjs").Unit): number | import("mathjs").Matrix | import("mathjs").MathArray;
    sparse(data?: import("mathjs").Matrix | import("mathjs").MathArray | undefined, dataType?: string | undefined): import("mathjs").Matrix;
    string(value: any): string | import("mathjs").Matrix | import("mathjs").MathArray;
    unit(unit: string): import("mathjs").Unit;
    unit(value: number, unit: string): import("mathjs").Unit;
    createUnit(name: string, definition?: string | import("mathjs").UnitDefinition | undefined, options?: import("mathjs").CreateUnitOptions | undefined): import("mathjs").Unit;
    createUnit(units: Record<string, string | import("mathjs").UnitDefinition>, options?: import("mathjs").CreateUnitOptions | undefined): import("mathjs").Unit;
    compile(expr: import("mathjs").MathExpression): import("mathjs").EvalFunction;
    compile(exprs: import("mathjs").MathExpression[]): import("mathjs").EvalFunction[];
    eval(expr: import("mathjs").MathExpression | import("mathjs").MathExpression[], scope?: any): any;
    help(search: any): import("mathjs").Help;
    parse(expr: import("mathjs").MathExpression, options?: any): import("mathjs").MathNode;
    parse(exprs: import("mathjs").MathExpression[], options?: any): import("mathjs").MathNode[];
    parser(): import("mathjs").Parser;
    distance(x: import("mathjs").MathType, y: import("mathjs").MathType): number | import("mathjs").BigNumber;
    intersect(w: import("mathjs").Matrix | import("mathjs").MathArray, x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").Matrix | import("mathjs").MathArray, z: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").MathArray;
    and(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray, y: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    not(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    or(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray, y: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    xor(x: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray, y: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    concat(...args: (number | import("mathjs").Matrix | import("mathjs").MathArray)[]): import("mathjs").Matrix | import("mathjs").MathArray;
    cross(x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix;
    det(x: import("mathjs").Matrix | import("mathjs").MathArray): number;
    diag(X: import("mathjs").Matrix | import("mathjs").MathArray, format?: string | undefined): import("mathjs").Matrix;
    diag(X: import("mathjs").Matrix | import("mathjs").MathArray, k: number | import("mathjs").BigNumber, format?: string | undefined): import("mathjs").Matrix;
    dot(x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").Matrix | import("mathjs").MathArray): number;
    eye(n: number | number[], format?: string | undefined): import("mathjs").Matrix;
    eye(m: number, n: number, format?: string | undefined): import("mathjs").Matrix;
    flatten(x: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    inv(x: number | import("mathjs").Matrix | import("mathjs").Complex | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").Complex | import("mathjs").MathArray;
    ones(n: number | number[], format?: string | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    ones(m: number, n: number, format?: string | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    range(str: string, includeEnd?: boolean | undefined): import("mathjs").Matrix;
    range(start: number | import("mathjs").BigNumber, end: number | import("mathjs").BigNumber, includeEnd?: boolean | undefined): import("mathjs").Matrix;
    range(start: number | import("mathjs").BigNumber, end: number | import("mathjs").BigNumber, step: number | import("mathjs").BigNumber, includeEnd?: boolean | undefined): import("mathjs").Matrix;
    resize(x: import("mathjs").Matrix | import("mathjs").MathArray, size: import("mathjs").Matrix | import("mathjs").MathArray, defaultValue?: string | number | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    size(x: string | number | boolean | import("mathjs").Matrix | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    squeeze(x: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    subset(value: string | import("mathjs").Matrix | import("mathjs").MathArray, index: import("mathjs").Index, replacement?: any, defaultValue?: any): string | import("mathjs").Matrix | import("mathjs").MathArray;
    trace(x: import("mathjs").Matrix | import("mathjs").MathArray): number;
    transpose(x: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    zeros(n: number | number[], format?: string | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    zeros(m: number, n: number, format?: string | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    distribution(name: string): import("mathjs").Distribution;
    factorial(n: number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").MathArray;
    gamma(n: number | import("mathjs").Matrix | import("mathjs").MathArray): number | import("mathjs").Matrix | import("mathjs").MathArray;
    kldivergence(x: import("mathjs").Matrix | import("mathjs").MathArray, y: import("mathjs").Matrix | import("mathjs").MathArray): number;
    multinomial(a: number[] | import("mathjs").BigNumber[]): number | import("mathjs").BigNumber;
    permutations(n: number | import("mathjs").BigNumber, k?: number | import("mathjs").BigNumber | undefined): number | import("mathjs").BigNumber;
    pickRandom(array: number[]): number;
    random(min?: number | undefined, max?: number | undefined): number;
    random(size: import("mathjs").Matrix | import("mathjs").MathArray, min?: number | undefined, max?: number | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    randomInt(min: number, max?: number | undefined): number;
    randomInt(size: import("mathjs").Matrix | import("mathjs").MathArray, min?: number | undefined, max?: number | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    compare(x: import("mathjs").MathType, y: import("mathjs").MathType): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").MathArray;
    deepEqual(x: import("mathjs").MathType, y: import("mathjs").MathType): number | import("mathjs").Matrix | import("mathjs").BigNumber | import("mathjs").Fraction | import("mathjs").Complex | import("mathjs").Unit | import("mathjs").MathArray;
    equal(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    larger(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    largerEq(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    smaller(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    smallerEq(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    unequal(x: import("mathjs").MathType, y: import("mathjs").MathType): boolean | import("mathjs").Matrix | import("mathjs").MathArray;
    max(...args: import("mathjs").MathType[]): any;
    max(A: import("mathjs").Matrix | import("mathjs").MathArray, dim?: number | undefined): any;
    mean(...args: import("mathjs").MathType[]): any;
    mean(A: import("mathjs").Matrix | import("mathjs").MathArray, dim?: number | undefined): any;
    median(...args: import("mathjs").MathType[]): any;
    min(...args: import("mathjs").MathType[]): any;
    min(A: import("mathjs").Matrix | import("mathjs").MathArray, dim?: number | undefined): any;
    mode(...args: import("mathjs").MathType[]): any;
    prod(...args: import("mathjs").MathType[]): any;
    quantileSeq(A: import("mathjs").Matrix | import("mathjs").MathArray, prob: number | import("mathjs").BigNumber | import("mathjs").MathArray, sorted?: boolean | undefined): number | import("mathjs").BigNumber | import("mathjs").Unit | import("mathjs").MathArray;
    std(array: import("mathjs").Matrix | import("mathjs").MathArray, normalization?: string | undefined): number;
    sum(...args: (number | import("mathjs").BigNumber | import("mathjs").Fraction)[]): any;
    sum(array: import("mathjs").Matrix | import("mathjs").MathArray): any;
    var(...args: (number | import("mathjs").BigNumber | import("mathjs").Fraction)[]): any;
    var(array: import("mathjs").Matrix | import("mathjs").MathArray, normalization?: string | undefined): any;
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
    atan2(y: import("mathjs").Matrix | import("mathjs").MathArray, x: import("mathjs").Matrix | import("mathjs").MathArray): import("mathjs").Matrix | import("mathjs").MathArray;
    atanh(x: number): number;
    atanh(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
    atanh(x: import("mathjs").MathArray): import("mathjs").MathArray;
    atanh(x: import("mathjs").Matrix): import("mathjs").Matrix;
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
    cos(x: number | import("mathjs").Unit): number;
    cos(x: import("mathjs").BigNumber): import("mathjs").BigNumber;
    cos(x: import("mathjs").Complex): import("mathjs").Complex;
    cos(x: import("mathjs").MathArray): import("mathjs").MathArray;
    cos(x: import("mathjs").Matrix): import("mathjs").Matrix;
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
    to(x: import("mathjs").Matrix | import("mathjs").Unit | import("mathjs").MathArray, unit: string | import("mathjs").Unit): import("mathjs").Matrix | import("mathjs").Unit | import("mathjs").MathArray;
    clone(x: any): any;
    filter(x: import("mathjs").Matrix | import("mathjs").MathArray, test: RegExp | ((item: any) => boolean)): import("mathjs").Matrix | import("mathjs").MathArray;
    forEach: (x: import("mathjs").Matrix | import("mathjs").MathArray, callback: (item: any) => any) => void;
    format(value: any, options?: number | import("mathjs").FormatOptions | ((item: any) => string) | undefined): string;
    isInteger(x: any): boolean;
    isNegative(x: any): boolean;
    isNumeric(x: any): boolean;
    isPositive(x: any): boolean;
    isZero(x: any): boolean;
    map(x: import("mathjs").Matrix | import("mathjs").MathArray, callback: (item: any) => any): import("mathjs").Matrix | import("mathjs").MathArray;
    partitionSelect(x: import("mathjs").Matrix | import("mathjs").MathArray, k: number, compare?: string | ((a: any, b: any) => number) | undefined): any;
    print: (template: string, values: any, precision?: number | undefined) => void;
    sort(x: import("mathjs").Matrix | import("mathjs").MathArray, compare?: string | ((a: any, b: any) => number) | undefined): import("mathjs").Matrix | import("mathjs").MathArray;
    typeof(x: any): string;
};
export default _default;
