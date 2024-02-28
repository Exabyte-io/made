export class InvalidLineError extends Error {
    constructor(num: any, content: any);
    lineNumber: any;
    content: any;
}
