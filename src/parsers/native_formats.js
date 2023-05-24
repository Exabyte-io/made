/**
 * @summary Function to handle coversion from native format
 * @returns {Object}
 */
function convertFromNative(format, text) {
    switch (format) {
        case "json":
            return JSON.parse(text);
        case "xyz":
            return {"message": "XYZ not supported yet"};
        case "poscar":
            return {"message": "POSCAR not supported yet"};
    }
}

export default convertFromNative;