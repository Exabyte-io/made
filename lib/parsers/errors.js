"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.InvalidLineError = void 0;
class InvalidLineError extends Error {
    constructor(num, content) {
        super(`Invalid line: ${num}`);
        this.lineNumber = num;
        this.content = content;
    }
}
exports.InvalidLineError = InvalidLineError;
//# sourceMappingURL=errors.js.map