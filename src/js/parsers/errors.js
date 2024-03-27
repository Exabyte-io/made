export class InvalidLineError extends Error {
    constructor(num, content) {
        super(`Invalid line: ${num}`);
        this.lineNumber = num;
        this.content = content;
    }
}
