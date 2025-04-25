export class InvalidLineError extends Error {
    constructor(num: number, content: string) {
        super(`Invalid line: ${num}`);
        console.log(`Invalid line: ${num}, content: ${content}`);
    }
}
