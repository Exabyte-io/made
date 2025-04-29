interface KPointStep {
    point: string;
    steps: number;
}
interface PathsType {
    [key: string]: KPointStep[];
}
export declare const paths: PathsType;
export {};
