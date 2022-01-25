export function getLineFromContent(fileContent, lineNumber) {
    return fileContent.split(/\r?\n/)[lineNumber];
}

export default {
    getLineFromContent,
};
