export function getLineFromContent(fileContent, lineNumber) {
    const lineString = fileContent.split(/\r?\n/)[lineNumber];
    return lineString.trim();
}
